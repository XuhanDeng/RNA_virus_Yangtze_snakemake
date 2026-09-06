#!/usr/bin/env python3
"""
Convert merged read_count TSV to FastSpar OTU table format.
All OTUs are kept regardless of prevalence.
"""

import argparse
import os
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Prepare read_count TSV for FastSpar."
    )
    parser.add_argument("--input",              required=True,
                        help="Merged *read_count.tsv file")
    parser.add_argument("--output",             required=True,
                        help="Output OTU table TSV for fastspar --otu_table")
    parser.add_argument("--min-samples-others", type=int, default=0,
                        help="kept for compatibility, no longer used")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    rc_cols = [c for c in df.columns if c.endswith("_read_count")]
    if not rc_cols:
        raise SystemExit("No *_read_count columns found in input.")

    label_col = df.columns[0]
    print(f"Label column : {label_col}  ({len(df)} OTUs, {len(rc_cols)} samples)",
          file=sys.stderr)

    # keep only label + count columns, rename
    out = df[[label_col] + rc_cols].copy()
    out = out.rename(columns={label_col: "OTU_id"})
    out = out.rename(columns={c: c.replace("_read_count", "") for c in rc_cols})

    sample_cols = [c.replace("_read_count", "") for c in rc_cols]
    for c in sample_cols:
        out[c] = pd.to_numeric(out[c], errors="coerce").fillna(0).astype(int)

    print(f"OTUs: {len(out)}", file=sys.stderr)
    print(f"Output: {len(out)} rows -> {args.output}", file=sys.stderr)
    out.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
