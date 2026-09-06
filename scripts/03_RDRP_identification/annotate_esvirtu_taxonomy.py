#!/usr/bin/env python3
"""
Annotate esvirtu-reference protein sequences with taxonomy looked up directly
from the ESvirtu metadata table (no DIAMOND search needed -- the taxonomy is
already known per accession in that table).

Protein headers embed the source accession as "esvirtu_{ACCESSION}", e.g.:
    >esvirtu_KX783432.1_frame=3_RdRp_245-2718 ...        (RdRpCATCH-style)
    >ORF2_esvirtu_AB162006.1:3558:5069 ...                (LucaProt-style)
The accession is the substring between "esvirtu_" and the next "_frame=" or
":" delimiter.

The metadata TSV (esvirtu_reference.metadata.tsv) has no header row; columns
used (1-indexed):
    1  accession
    5  k__Kingdom
    6  p__Phylum
    7  c__Class
    8  o__Order
    9  f__Family
    10 g__Genus
    11 s__Species

Output TSV matches the column shape expected by split_fasta_by_phylum.py
(qseqid, hit_rank, sseqid, Phylum, Class, Order, Family), with hit_rank
always "1" and sseqid always blank since this is a direct lookup, not a
DIAMOND search result.
"""

import argparse
import csv
import re
import sys

ACC_RE = re.compile(r"esvirtu_(.+?)(?:_frame=|:\d)")

TAX_COLS = {
    "Phylum": 6,
    "Class":  7,
    "Order":  8,
    "Family": 9,
    "Genus":  10,
    "Species": 11,
}

OUT_COLS = ["qseqid", "hit_rank", "sseqid", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def strip_prefix(value: str) -> str:
    """Strip a GTDB-style 'x__' prefix, e.g. 'p__Pisuviricota' -> 'Pisuviricota'."""
    value = value.strip()
    return re.sub(r"^[a-z]__", "", value)


def load_metadata(path):
    """accession -> {rank_name: value}."""
    tax = {}
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(TAX_COLS.values()) - 1:
                continue
            accession = cols[0].strip()
            tax[accession] = {
                rank: strip_prefix(cols[idx - 1])
                for rank, idx in TAX_COLS.items()
            }
    return tax


def iter_fasta_headers(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].rstrip("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta",    required=True,
                         help="Protein FASTA to annotate (e.g. combined_full.faa)")
    parser.add_argument("--metadata", required=True,
                         help="esvirtu_reference.metadata.tsv")
    parser.add_argument("--output",   required=True)
    args = parser.parse_args()

    print("Loading ESvirtu metadata taxonomy...", file=sys.stderr)
    tax = load_metadata(args.metadata)
    print(f"  {len(tax)} accessions loaded", file=sys.stderr)

    n_total = n_matched = n_missing = 0
    with open(args.output, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=OUT_COLS, delimiter="\t")
        writer.writeheader()

        for header in iter_fasta_headers(args.fasta):
            n_total += 1
            seq_id = header.split()[0]
            m = ACC_RE.search(seq_id)
            accession = m.group(1) if m else None

            row = {col: "" for col in OUT_COLS}
            row["qseqid"] = seq_id
            row["hit_rank"] = "1"
            row["sseqid"] = ""

            entry = tax.get(accession) if accession else None
            if entry is None:
                n_missing += 1
                print(f"  MISSING taxonomy: {seq_id} (accession={accession})",
                      file=sys.stderr)
            else:
                n_matched += 1
                for rank in ("Phylum", "Class", "Order", "Family", "Genus", "Species"):
                    row[rank] = entry.get(rank, "")

            writer.writerow(row)

    print(f"Annotated {n_total} sequences: {n_matched} matched, "
          f"{n_missing} missing taxonomy", file=sys.stderr)


if __name__ == "__main__":
    main()
