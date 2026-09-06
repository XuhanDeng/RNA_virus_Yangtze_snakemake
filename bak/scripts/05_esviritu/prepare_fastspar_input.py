#!/usr/bin/env python3
"""
Convert merged read_count TSV to FastSpar OTU table format.

OTUs present in fewer than --min-samples-others samples are collapsed
into a single 'others' row (read counts summed across samples).
OTUs meeting the threshold are kept as individual rows.
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
    parser.add_argument("--min-samples-others", type=int, default=3,
                        help="OTUs present in fewer than N samples are collapsed "
                             "into an 'others' row (default: 3)")
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

    # split into prevalent (keep) and rare (collapse to others)
    presence = (out[sample_cols] > 0).sum(axis=1)
    keep_mask = presence >= args.min_samples_others
    rare_mask  = ~keep_mask

    n_keep = keep_mask.sum()
    n_rare = rare_mask.sum()
    print(f"OTUs kept as individual rows : {n_keep}", file=sys.stderr)
    print(f"OTUs collapsed into 'others' : {n_rare}", file=sys.stderr)

    kept = out[keep_mask].copy()

    if rare_mask.any():
        others_counts = out.loc[rare_mask, sample_cols].sum(axis=0).astype(int)
        others_row = pd.DataFrame([["others"] + others_counts.tolist()],
                                  columns=["OTU_id"] + sample_cols)
        result = pd.concat([kept, others_row], ignore_index=True)
    else:
        result = kept

    print(f"Output: {len(result)} rows (including 'others') -> {args.output}",
          file=sys.stderr)
    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
