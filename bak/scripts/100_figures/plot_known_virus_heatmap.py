#!/usr/bin/env python3
"""
Heatmap of per-sample abundance for known ESvirtu-reference viruses (seq_id
starting with "t__"), one row per virus (no grouping).

Reads a recalculated RPKMF or TPM table (rna_virus_recalc_rpkmf.tsv /
rna_virus_recalc_tpm.tsv), keeps only rows where seq_id starts with "t__",
labels each row "<seq_id> (<merged_host>)" (or just "<seq_id>" if
merged_host is empty), ranks viruses by total abundance across all samples,
and draws a virus x sample heatmap (log1p-scaled color) for the top-N
viruses. Samples ordered left-to-right by trailing numeric sample index.
"""

import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")
# keep SVG text as real <text> elements (editable in Illustrator),
# instead of being converted to paths.
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

KNOWN_PREFIX = "t__"


def sample_sort_key(value_col: str, suffix: str):
    """Sort by the trailing numeric sample index (e.g. ..._42_rpkmf -> 42)."""
    match = re.search(rf"(\d+){re.escape(suffix)}$", value_col)
    return int(match.group(1)) if match else float("inf")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True,
                         help="rna_virus_recalc_rpkmf.tsv or rna_virus_recalc_tpm.tsv")
    parser.add_argument("--value-suffix", required=True,
                         help="suffix identifying per-sample value columns (e.g. _rpkmf or _tpm)")
    parser.add_argument("--top-n", type=int, required=True,
                         help="number of top known viruses (by total abundance) to show in the heatmap")
    parser.add_argument("--output", required=True, help="output figure path (.svg/.png)")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    value_cols = [c for c in df.columns if c.endswith(args.value_suffix)]
    if not value_cols:
        sys.exit(f"No {args.value_suffix} columns found in input.")
    value_cols = sorted(value_cols, key=lambda c: sample_sort_key(c, args.value_suffix))

    seq_id = df.get("seq_id")
    merged_host = df.get("merged_host")
    if seq_id is None:
        sys.exit("No seq_id column found in input.")
    if merged_host is None:
        sys.exit("No merged_host column found in input.")

    seq_id = seq_id.astype(str)
    df = df[seq_id.str.startswith(KNOWN_PREFIX)].copy()
    if df.empty:
        sys.exit(f"No rows with seq_id starting with '{KNOWN_PREFIX}' found in input.")

    merged_host = df["merged_host"].fillna("").astype(str).str.strip()
    df["label"] = df["seq_id"].where(
        merged_host == "", df["seq_id"] + " (" + merged_host + ")"
    )

    for c in value_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    named = df.set_index("label")[value_cols].copy()
    samples = [c[: -len(args.value_suffix)] for c in value_cols]
    named.columns = samples

    # rank viruses by total abundance across all samples, descending; keep top-N
    totals = named.sum(axis=1).sort_values(ascending=False)
    top_viruses = totals.head(args.top_n).index.tolist()
    heat_named = named.loc[top_viruses]

    # ── plot ──────────────────────────────────────────────────────────────────
    data = heat_named.to_numpy(dtype=float)
    log_data = np.log1p(data)

    n_rows, n_cols = data.shape
    fig, ax = plt.subplots(figsize=(max(10, n_cols * 0.35), max(4, n_rows * 0.4 + 2)))

    im = ax.imshow(log_data, aspect="auto", cmap="viridis")

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(samples, rotation=90)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(heat_named.index, fontsize=7)

    ax.set_xlabel("Sample")
    ax.set_ylabel("known virus (seq_id (merged_host))")
    ax.set_title(f"Top {args.top_n} known-virus abundance per sample (log1p scale)")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("log1p(abundance)")

    fig.tight_layout()
    fig.savefig(args.output)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
