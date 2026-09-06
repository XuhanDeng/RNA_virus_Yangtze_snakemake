#!/usr/bin/env python3
"""
Filter final RdRp protein candidates by Diamond BLASTp hit against nr,
restricted to the Riboviria taxon.

The calling rule already ran:
    diamond blastp --taxonlist <riboviria_taxid> --evalue <threshold> -k 1 ...
so `--diamond` only ever contains, per query, its single best hit *within
Riboviria* (if any) passing the e-value threshold. A query with no row in
that file therefore has no significant RNA virus match in nr at all, and is
dropped -- this catches endogenous viral elements and other false positives
that pass the motif filter but aren't of confirmed viral origin.

Reads:
  --diamond   diamond blastp --outfmt 6 qseqid sseqid evalue bitscore staxids
              output, already restricted to Riboviria hits (--taxonlist)
  --fasta     the protein FASTA that was searched (e.g. RdRpCATCH.faa or
              LucaProt.faa) -- every sequence in here is a query

Writes:
  --output-fasta  same records as --fasta, restricted to queries with a
                  Riboviria hit
  --output-table  one row per query: kept/dropped, top hit accession,
                  taxids, e-value, reason
"""

import argparse
import os
import sys

import pandas as pd

DIAMOND_COLS = ["qseqid", "sseqid", "evalue", "bitscore", "staxids"]


def _load_best_hits(diamond_path: str) -> pd.DataFrame:
    if not os.path.exists(diamond_path) or os.path.getsize(diamond_path) == 0:
        return pd.DataFrame(columns=DIAMOND_COLS)
    df = pd.read_csv(diamond_path, sep="\t", header=None, names=DIAMOND_COLS, dtype=str)
    if df.empty:
        return df
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    # -k 1 already restricts diamond's own output to one hit per query, but
    # sort/dedup defensively in case that ever changes.
    df = df.sort_values(["qseqid", "evalue"]).drop_duplicates("qseqid", keep="first")
    return df


def _read_fasta_ids(fasta_path: str):
    ids = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def _write_filtered_fasta(fasta_path: str, keep_ids: set, output_path: str):
    keep = False
    with open(fasta_path) as fin, open(output_path, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                qid = line[1:].strip().split()[0]
                keep = qid in keep_ids
            if keep:
                fout.write(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diamond", required=True,
                         help="diamond blastp --outfmt 6 qseqid sseqid evalue bitscore staxids, "
                              "already restricted to Riboviria hits via --taxonlist")
    parser.add_argument("--fasta", required=True,
                         help="protein FASTA that was searched against nr")
    parser.add_argument("--output-fasta", required=True)
    parser.add_argument("--output-table", required=True)
    args = parser.parse_args()

    all_ids = _read_fasta_ids(args.fasta)
    print(f"Query proteins in {args.fasta}: {len(all_ids)}", file=sys.stderr)

    best_hits = _load_best_hits(args.diamond)
    hit_by_query = {row["qseqid"]: row for _, row in best_hits.iterrows()}

    rows = []
    keep_ids = set()
    for qid in all_ids:
        hit = hit_by_query.get(qid)
        if hit is None:
            rows.append({"qseqid": qid, "kept": False, "sseqid": None,
                         "staxids": None, "evalue": None,
                         "reason": "no_riboviria_hit"})
            continue

        rows.append({
            "qseqid":  qid,
            "kept":    True,
            "sseqid":  hit.get("sseqid"),
            "staxids": hit.get("staxids"),
            "evalue":  hit.get("evalue"),
            "reason":  "riboviria_hit",
        })
        keep_ids.add(qid)

    out = pd.DataFrame(rows)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_table)), exist_ok=True)
    out.to_csv(args.output_table, sep="\t", index=False)

    os.makedirs(os.path.dirname(os.path.abspath(args.output_fasta)), exist_ok=True)
    _write_filtered_fasta(args.fasta, keep_ids, args.output_fasta)

    print(f"Kept {len(keep_ids)} / {len(all_ids)} sequences (Riboviria hit in nr)", file=sys.stderr)
    print(f"  -> {args.output_fasta}", file=sys.stderr)
    print(f"  -> {args.output_table}", file=sys.stderr)


if __name__ == "__main__":
    main()
