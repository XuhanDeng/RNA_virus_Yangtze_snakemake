#!/usr/bin/env python3
"""
Heatmap of per-sample abundance grouped by host (minimum_host_lineage,
falling back to merged_host when the lineage is empty).

Reads a recalculated RPKMF or TPM table (rna_virus_recalc_rpkmf.tsv /
rna_virus_recalc_tpm.tsv), sums each sample's value columns grouped by a
per-row host label: minimum_host_lineage if present, else merged_host,
else the row is dropped (true unknown). Writes --output-table: the full
abundance table for every clustered host (no cutoff), sorted by descending
total abundance. The heatmap figure itself only shows the top-N hosts
(--top-n) by total abundance. Samples ordered left-to-right by trailing
numeric sample index.
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
                         help="number of top hosts (by total abundance) to show in the heatmap")
    parser.add_argument("--output", required=True, help="output figure path (.svg/.png)")
    parser.add_argument("--output-table", required=True,
                         help="output CSV path with per-sample, per-host abundance (ALL hosts, no cutoff)")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", low_memory=False)

    value_cols = [c for c in df.columns if c.endswith(args.value_suffix)]
    if not value_cols:
        sys.exit(f"No {args.value_suffix} columns found in input.")
    value_cols = sorted(value_cols, key=lambda c: sample_sort_key(c, args.value_suffix))

    lineage = df.get("minimum_host_lineage")
    merged_host = df.get("merged_host")
    if lineage is None:
        sys.exit("No minimum_host_lineage column found in input.")
    if merged_host is None:
        sys.exit("No merged_host column found in input.")

    lineage = lineage.fillna("").astype(str).str.strip()
    merged_host = merged_host.fillna("").astype(str).str.strip()

    # host label: minimum_host_lineage if present, else merged_host, else drop (true unknown)
    df["host"] = lineage.where(lineage != "", merged_host)
    df = df[df["host"] != ""]

    for c in value_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    # sum abundance per (sample, host)
    grouped = df.groupby("host")[value_cols].sum()

    samples = [c[: -len(args.value_suffix)] for c in value_cols]
    grouped_named = grouped.copy()
    grouped_named.columns = samples

    # rank hosts by total abundance across all samples, descending
    totals_per_host = grouped_named.sum(axis=1).sort_values(ascending=False)
    grouped_named = grouped_named.loc[totals_per_host.index]

    # ── CSV table: full, every clustered host, no cutoff, wide format ───────────
    # one row per host, one column per sample, sorted by descending total abundance
    full_table = grouped_named.reset_index()
    full_table.to_csv(args.output_table, index=False)
    print(f"Wrote {args.output_table}", file=sys.stderr)

    # heatmap figure only shows the top-N hosts by total abundance
    top_hosts = totals_per_host.head(args.top_n).index.tolist()
    heat_named = grouped_named.loc[top_hosts]

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
    ax.set_ylabel("host (minimum_host_lineage / merged_host)")
    ax.set_title(f"Top {args.top_n} host abundance per sample (log1p scale)")

    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("log1p(abundance)")

    fig.tight_layout()
    fig.savefig(args.output)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
