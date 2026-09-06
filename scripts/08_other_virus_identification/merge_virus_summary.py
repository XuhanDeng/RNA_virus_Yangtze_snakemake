#!/usr/bin/env python3
"""Merge per-sample GeNomad virus_summary TSVs, adding a 'sample' column."""

import argparse
import pathlib
import sys

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("inputs", nargs="+", help="Per-sample virus_summary.tsv files")
    args = parser.parse_args()

    dfs = []
    for f in args.inputs:
        sample = pathlib.Path(f).name.replace("_virus_filter_summary.tsv", "")
        df = pd.read_csv(f, sep="\t", dtype=str)
        df.insert(0, "sample", sample)
        dfs.append(df)

    merged = pd.concat(dfs, ignore_index=True)
    merged.to_csv(args.output, sep="\t", index=False)
    print(f"Merged {len(merged)} rows from {len(dfs)} samples -> {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
