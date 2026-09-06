#!/usr/bin/env python3
"""
Annotate known_unknown_pair TSV files with contig_tag.

Reads rna_virus_contig_table.tsv and maps each contig to its contig_tag:
  RdRpCATCH, LucaProt,
  bait_RdRpCATCH, bait_LucaProt, bait_RVMT,
  unknown — contig not found in rna_virus_contig_table
"""

import argparse
import os
import sys

import pandas as pd


def _load_tag_map(contig_table: str) -> dict:
    df = pd.read_csv(contig_table, sep="\t", dtype=str,
                     usecols=lambda c: c in {"contig_id", "contig_tag"})
    df = df.dropna(subset=["contig_id", "contig_tag"]).drop_duplicates("contig_id")

    counts = df["contig_tag"].value_counts().to_dict()
    print(f"contig_tag counts from {contig_table}: {counts}", file=sys.stderr)

    return dict(zip(df["contig_id"], df["contig_tag"]))


def main():
    parser = argparse.ArgumentParser(
        description="Annotate known_unknown_pair TSV with contig_tag from rna_virus_contig_table."
    )
    parser.add_argument("--input",         required=True, help="Input TSV (known_unknown_pair)")
    parser.add_argument("--output",        required=True, help="Output TSV")
    parser.add_argument("--contig-table",  required=True, help="Path to rna_virus_contig_table.tsv")
    args = parser.parse_args()

    contig_to_tag = _load_tag_map(args.contig_table)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    if df.empty or "contig_name" not in df.columns:
        df["contig_tag"] = pd.Series(dtype=str)
        df.to_csv(args.output, sep="\t", index=False)
        print(f"Input empty or missing contig_name — written empty output: {args.output}",
              file=sys.stderr)
        return

    df["contig_tag"] = df["contig_name"].map(contig_to_tag).fillna("unknown")

    counts = df["contig_tag"].value_counts().to_dict()
    print(f"Annotated {len(df)} rows: {counts} -> {args.output}", file=sys.stderr)

    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
