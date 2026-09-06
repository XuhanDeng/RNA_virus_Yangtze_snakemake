#!/usr/bin/env python3
"""
Concatenate all per-sample RdRp merged TSVs into one final table and write
contig list files.

Each per-sample TSV (from merge_rdrp_results.py) has one row per contig with
a `source` column (RdRpCATCH or LucaProt) indicating which tool's motif-
confirmed region was kept for that contig.

List files written
-------------------
  RdRpCATCH.txt — contigs whose kept region came from RdRpCATCH
  LucaProt.txt  — contigs whose kept region came from LucaProt
  any_rdrp.txt  — union of both (== all contigs in the combined table)
"""

import argparse
import os
import sys

import pandas as pd


def _write_list(path: str, contigs: pd.Series):
    contigs = contigs.dropna().drop_duplicates()
    contigs.to_csv(path, index=False, header=False)
    print(f"  {os.path.basename(path)}: {len(contigs)} contigs", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Combine per-sample RdRp merged TSVs and write contig list files."
    )
    parser.add_argument("--inputs", nargs="+", required=True,
                        help="Per-sample *_rdrp_merged.tsv files")
    parser.add_argument("--input-regions", nargs="+", required=False, default=None,
                        help="Per-sample *_rdrp_regions.tsv files (pre-longest-only candidates)")
    parser.add_argument("--output", required=True,
                        help="Output combined contig TSV")
    parser.add_argument("--output-regions", required=False, default=None,
                        help="Output combined region TSV")
    parser.add_argument("--list-dir", required=True,
                        help="Directory to write contig list files")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    os.makedirs(args.list_dir, exist_ok=True)

    parts = []
    for path in args.inputs:
        sample = os.path.basename(path).replace("_rdrp_merged.tsv", "")
        df = pd.read_csv(path, sep="\t", dtype=str)
        df.insert(0, "sample", sample)
        parts.append(df)
        print(f"  {sample}: {len(df)} contigs", file=sys.stderr)

    combined = pd.concat(parts, ignore_index=True)
    combined.to_csv(args.output, sep="\t", index=False)
    print(f"Combined: {len(combined)} total rows from {len(parts)} samples -> {args.output}",
          file=sys.stderr)

    if args.input_regions and args.output_regions:
        region_parts = []
        for path in args.input_regions:
            sample = os.path.basename(path).replace("_rdrp_regions.tsv", "")
            df = pd.read_csv(path, sep="\t", dtype=str)
            df.insert(0, "sample", sample)
            region_parts.append(df)
        combined_regions = pd.concat(region_parts, ignore_index=True)
        os.makedirs(os.path.dirname(os.path.abspath(args.output_regions)), exist_ok=True)
        combined_regions.to_csv(args.output_regions, sep="\t", index=False)
        print(f"Combined regions: {len(combined_regions)} rows -> {args.output_regions}",
              file=sys.stderr)

    is_rc = combined["source"] == "RdRpCATCH"
    is_lp = combined["source"] == "LucaProt"

    print("Writing list files:", file=sys.stderr)
    _write_list(os.path.join(args.list_dir, "RdRpCATCH.txt"), combined.loc[is_rc, "contig"])
    _write_list(os.path.join(args.list_dir, "LucaProt.txt"),  combined.loc[is_lp, "contig"])
    _write_list(os.path.join(args.list_dir, "any_rdrp.txt"),  combined["contig"])


if __name__ == "__main__":
    main()
