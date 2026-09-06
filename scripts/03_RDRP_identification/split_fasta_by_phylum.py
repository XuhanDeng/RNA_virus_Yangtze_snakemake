#!/usr/bin/env python3
"""
Split a protein FASTA into per-taxon FASTA files using the top DIAMOND-vs-RVMT
hit (hit_rank == 1) taxonomy call for each sequence, at a chosen rank
(Phylum, Class, Order, or Family).

Sequences with no hit_rank==1 row (no significant DIAMOND hit at all) or a
blank value at the chosen rank are written to "Unclassified.faa".

If --rvmt-fasta is given, each query's top-hit RVMT reference sequence
(hit_rank==1 sseqid) is also extracted and written to a separate
"{taxon}_ref.faa" file, deduplicated (multiple queries can share the same
top RVMT hit).

Usage:
    python split_fasta_by_phylum.py \
        --fasta      combined_full_c90.faa \
        --annotated  full_length_diamond_rvmt_annotated.tsv \
        --outdir     22_phylum_cluster/full_length \
        --rank       Phylum \
        --rvmt-fasta RdRPs_ali822x.faa
"""

import argparse
import csv
import os
import re
import sys

RANKS = ["Phylum", "Class", "Order", "Family"]


def read_fasta(path):
    """Yield (header, seq) tuples. header is everything after '>' on its line."""
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            elif header is not None:
                seq_lines.append(line)
        if header is not None:
            yield header, "".join(seq_lines)


def wrap(seq, width=60):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def safe_filename(taxon):
    if not taxon:
        return "Unclassified"
    return re.sub(r"[^A-Za-z0-9._-]", "_", taxon)


def load_top_hits(annotated_tsv, rank):
    """qseqid -> (taxon_at_rank, sseqid), using only hit_rank == 1 rows."""
    top_hits = {}
    with open(annotated_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row.get("hit_rank") != "1":
                continue
            taxon = row.get(rank, "").strip()
            sseqid = row.get("sseqid", "").strip()
            top_hits[row["qseqid"]] = (taxon, sseqid)
    return top_hits


def index_fasta(path):
    """seq_id (first whitespace-delimited token of header) -> (header, seq)."""
    index = {}
    for header, seq in read_fasta(path):
        seq_id = header.split()[0]
        index[seq_id] = (header, seq)
    return index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta",      required=True)
    parser.add_argument("--annotated",  required=True)
    parser.add_argument("--outdir",     required=True)
    parser.add_argument("--rank",       default="Phylum", choices=RANKS,
                         help="Taxonomic rank column to split by (default: Phylum)")
    parser.add_argument("--rvmt-fasta", default=None,
                         help="If given, also extract each query's top-hit "
                              "RVMT reference sequence into {taxon}_ref.faa")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading top-hit {args.rank}/sseqid calls...", file=sys.stderr)
    top_hits = load_top_hits(args.annotated, args.rank)
    print(f"  {len(top_hits)} queries with a hit_rank==1 call.", file=sys.stderr)

    rvmt_index = None
    if args.rvmt_fasta:
        print("Indexing RVMT reference FASTA...", file=sys.stderr)
        rvmt_index = index_fasta(args.rvmt_fasta)
        print(f"  {len(rvmt_index)} RVMT sequences indexed.", file=sys.stderr)

    handles = {}
    counts = {}
    n_seqs = 0
    for header, seq in read_fasta(args.fasta):
        n_seqs += 1
        seq_id = header.split()[0]
        taxon, sseqid = top_hits.get(seq_id, ("", ""))
        taxon = taxon or "Unclassified"
        fname = safe_filename(taxon)

        if fname not in handles:
            handles[fname] = open(os.path.join(args.outdir, f"{fname}.faa"), "w")
            counts[fname] = 0

        handles[fname].write(f">{header}\n{wrap(seq)}\n")
        counts[fname] += 1

    for fh in handles.values():
        fh.close()

    print(f"Wrote {n_seqs} sequences across {len(handles)} {args.rank} files:", file=sys.stderr)
    for fname, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {fname}.faa: {n}", file=sys.stderr)

    if rvmt_index is None:
        return

    print(f"Extracting RVMT reference sequences per {args.rank}...", file=sys.stderr)
    ref_handles = {}
    ref_written = {}
    n_ref_missing = 0
    for header, seq in read_fasta(args.fasta):
        seq_id = header.split()[0]
        taxon, sseqid = top_hits.get(seq_id, ("", ""))
        if not sseqid:
            continue
        fname = safe_filename(taxon or "Unclassified")

        ref_written.setdefault(fname, set())
        if sseqid in ref_written[fname]:
            continue  # dedupe: multiple queries can share the same top RVMT hit

        rvmt_entry = rvmt_index.get(sseqid)
        if rvmt_entry is None:
            n_ref_missing += 1
            continue
        ref_header, ref_seq = rvmt_entry

        if fname not in ref_handles:
            ref_handles[fname] = open(os.path.join(args.outdir, f"{fname}_ref.faa"), "w")
        ref_handles[fname].write(f">{ref_header}\n{wrap(ref_seq)}\n")
        ref_written[fname].add(sseqid)

    for fh in ref_handles.values():
        fh.close()

    print(f"Wrote reference sequences across {len(ref_handles)} {args.rank} files:", file=sys.stderr)
    for fname, ids in sorted(ref_written.items(), key=lambda x: -len(x[1])):
        print(f"  {fname}_ref.faa: {len(ids)}", file=sys.stderr)
    if n_ref_missing:
        print(f"  Warning: {n_ref_missing} top-hit sseqids not found in RVMT fasta.", file=sys.stderr)


if __name__ == "__main__":
    main()
