#!/usr/bin/env python3
"""
Convert FastSpar wide-format correlation + p-value matrices to a
filtered long-format TSV matching the pipeline's standard output.

FastSpar output matrices:
  correlation.tsv  — square OTU x OTU matrix, values in [-1, 1]
  pvalues.tsv      — square OTU x OTU matrix, values in [0, 1]

Output columns: seed, target, r, p, r_group
Only pairs passing --r-threshold and --p-threshold are written.
Self-pairs (seed == target) are always excluded.
Each pair is reported once (seed < target alphabetically).
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd


def _parse_thresholds(text: str):
    vals = [float(p.strip()) for p in text.split(",") if p.strip()]
    if not vals:
        raise SystemExit("No valid thresholds provided.")
    return sorted(set(vals))


def main():
    parser = argparse.ArgumentParser(
        description="Convert FastSpar matrices to filtered long-format TSV."
    )
    parser.add_argument("--correlation", required=True,
                        help="FastSpar correlation matrix TSV")
    parser.add_argument("--pvalues",     required=True,
                        help="FastSpar p-value matrix TSV")
    parser.add_argument("--output",      required=True,
                        help="Output filtered long-format TSV")
    parser.add_argument("--thresholds",  default="0.6,0.7,0.8,0.9",
                        help="Comma-separated r thresholds for r_group column")
    parser.add_argument("--r-threshold", type=float, default=0.6,
                        help="Minimum |r| to keep a pair (default: 0.6)")
    parser.add_argument("--p-threshold", type=float, default=0.05,
                        help="Maximum p-value to keep a pair (default: 0.05)")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    thresholds = _parse_thresholds(args.thresholds)
    min_thr    = args.r_threshold

    corr = pd.read_csv(args.correlation, sep="\t", index_col=0)
    pval = pd.read_csv(args.pvalues,     sep="\t", index_col=0)

    otus = corr.index.tolist()
    n    = len(otus)
    print(f"Matrix size: {n} x {n}", file=sys.stderr)

    r_mat = corr.to_numpy(dtype=float)
    p_mat = pval.to_numpy(dtype=float)
    otus_arr = np.array(otus)

    # upper triangle only (i < j) to report each pair once, excluding diagonal
    i_idx, j_idx = np.triu_indices(n, k=1)
    r_vals = r_mat[i_idx, j_idx]
    p_vals = p_mat[i_idx, j_idx]

    mask = (r_vals >= min_thr) & (p_vals <= args.p_threshold)
    r_f  = r_vals[mask]
    p_f  = p_vals[mask]
    ia   = i_idx[mask]
    ib   = j_idx[mask]

    print(f"Pairs passing filters: {mask.sum()}", file=sys.stderr)

    if r_f.size == 0:
        pd.DataFrame(columns=["seed", "target", "r", "p", "r_group"]).to_csv(
            args.output, sep="\t", index=False)
        print(f"Written (empty): {args.output}", file=sys.stderr)
        return

    r_group = np.array([max(t for t in thresholds if rv >= t) for rv in r_f])

    out = pd.DataFrame({
        "seed":    otus_arr[ia],
        "target":  otus_arr[ib],
        "r":       r_f,
        "p":       p_f,
        "r_group": r_group,
    })

    out.to_csv(args.output, sep="\t", index=False)
    print(f"Written: {len(out)} pairs -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
