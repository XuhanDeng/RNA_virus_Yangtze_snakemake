#!/usr/bin/env python3
"""
Write per-source contig list .txt files from the ICTV merged TSV.

Reads the `source` column (RdRpCATCH or LucaProt) and writes RdRpCATCH.txt,
LucaProt.txt, plus any_rdrp.txt (union of both). Used by the ICTV pipeline
only; RNA samples use combine_rdrp_all_samples.py instead.
"""

import argparse
import os
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged", required=True,
                        help="ICTV_rdrp_merged.tsv from merge_rdrp_results.py")
    parser.add_argument("--outdir", required=True,
                        help="Directory to write list .txt files")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.merged, sep="\t", dtype=str)
    print(f"Loaded {len(df)} contigs from {args.merged}", file=sys.stderr)

    is_rc = df["source"] == "RdRpCATCH"
    is_lp = df["source"] == "LucaProt"

    for name, mask in [("RdRpCATCH", is_rc), ("LucaProt", is_lp)]:
        contigs = df.loc[mask, "contig"].dropna().drop_duplicates()
        path = os.path.join(args.outdir, f"{name}.txt")
        contigs.to_csv(path, index=False, header=False)
        print(f"  {name}: {len(contigs)} -> {path}", file=sys.stderr)

    any_path = os.path.join(args.outdir, "any_rdrp.txt")
    any_contigs = df["contig"].dropna().drop_duplicates()
    any_contigs.to_csv(any_path, index=False, header=False)
    print(f"  any_rdrp: {len(any_contigs)} -> {any_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
