#!/usr/bin/env python3
"""
Annotate known_unknown_pair TSV files with RdRP tier.

Reads tier_summary.tsv (output of motif_merge_filter.py) and maps each
contig to its tier:

  tier1a  — all 4 motifs, canonical ABCD order
  tier1b  — all 4 motifs, CABD (depermuted)
  tier2   — 3 motifs detected
  tier3   — identified by rdrpcatch/lucaprot but fewer than 3 motifs
  not_in_list — contig not found in tier_summary
"""

import argparse
import os
import sys

import pandas as pd


def _load_tier_map(tier_tsv: str) -> dict:
    df = pd.read_csv(tier_tsv, sep="\t", dtype=str,
                     usecols=lambda c: c in {"contig", "tier"})
    df = df.dropna(subset=["contig", "tier"]).drop_duplicates("contig")

    counts = df["tier"].value_counts().to_dict()
    print(f"Tier counts from {tier_tsv}: {counts}", file=sys.stderr)

    return dict(zip(df["contig"], df["tier"]))


def main():
    parser = argparse.ArgumentParser(
        description="Annotate known_unknown_pair TSV with rdrp_category (tier1a/1b/2/3)."
    )
    parser.add_argument("--input",    required=True, help="Input TSV (known_unknown_pair)")
    parser.add_argument("--output",   required=True, help="Output TSV")
    parser.add_argument("--tier-tsv", required=True, help="Path to tier_summary.tsv")
    args = parser.parse_args()

    contig_to_tier = _load_tier_map(args.tier_tsv)

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    if df.empty or "contig_name" not in df.columns:
        df["rdrp_category"] = pd.Series(dtype=str)
        df.to_csv(args.output, sep="\t", index=False)
        print(f"Input empty or missing contig_name — written empty output: {args.output}",
              file=sys.stderr)
        return

    df["rdrp_category"] = df["contig_name"].map(contig_to_tier).fillna("not_in_list")

    counts = df["rdrp_category"].value_counts().to_dict()
    print(f"Annotated {len(df)} rows: {counts} -> {args.output}", file=sys.stderr)

    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
