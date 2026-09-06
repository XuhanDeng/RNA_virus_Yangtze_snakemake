#!/usr/bin/env python3
"""
Bar plot of Homo sapiens contamination (%) per sample from Kraken2 reports.

For each sample's *_kraken2.report, extracts the percentage assigned to
Homo sapiens (S-rank, taxid 9606) and draws a sorted horizontal bar chart.
"""

import argparse
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd


def parse_kraken_report(path: str) -> float:
    """Return Homo sapiens percentage from a Kraken2 report, or 0.0 if absent."""
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            rank = parts[3].strip()
            name = parts[5].strip()
            if rank == "S" and name == "Homo sapiens":
                return float(parts[0].strip())
    return 0.0


def collect_data(kraken_dir: str) -> pd.DataFrame:
    """Scan kraken_dir for *_kraken2.report files and collect Homo sapiens %."""
    records = []
    pattern = os.path.join(kraken_dir, "*", "*_kraken2.report")
    paths = sorted(glob.glob(pattern))
    if not paths:
        sys.exit(f"ERROR: no *_kraken2.report files found under {kraken_dir}")
    for path in paths:
        sample = os.path.basename(os.path.dirname(path))
        pct = parse_kraken_report(path)
        records.append({"sample": sample, "human_pct": pct})
    return pd.DataFrame(records)


def plot(df: pd.DataFrame, outpath: str, threshold: float):
    # Sort by trailing numeric suffix (_1, _2, ... _48)
    import re as _re
    df["_sort_key"] = df["sample"].apply(
        lambda s: int(_re.search(r"_(\d+)$", s).group(1)) if _re.search(r"_(\d+)$", s) else 0
    )
    df = df.sort_values("_sort_key").drop(columns="_sort_key").reset_index(drop=True)

    n = len(df)
    fig_width = max(6, n * 0.28)
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    colors = ["#E6194B" if v >= threshold else "#4363D8" for v in df["human_pct"]]
    bars = ax.bar(df["sample"], df["human_pct"], color=colors, width=0.7, edgecolor="none")

    # Value labels above bars
    for bar, val in zip(bars, df["human_pct"]):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{val:.2f}%", va="bottom", ha="center", fontsize=6, rotation=90)

    # Threshold line
    ax.axhline(threshold, color="#E6194B", linestyle="--", linewidth=1, alpha=0.7)

    ax.set_ylabel("Homo sapiens reads (%)", fontsize=10)
    ax.set_title("Human Contamination per Sample (Kraken2)", fontsize=11, fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.1f%%"))
    plt.xticks(rotation=45, ha="right", fontsize=7)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#E6194B", label=f"≥ {threshold}% (high)"),
        Patch(facecolor="#4363D8", label=f"< {threshold}% (low)"),
        plt.Line2D([0], [0], color="#E6194B", linestyle="--", linewidth=1,
                   label=f"Threshold ({threshold}%)"),
    ]
    ax.legend(handles=legend_elements, fontsize=8, loc="upper right")

    plt.tight_layout()
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    print(f"Written → {outpath}", file=sys.stderr)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kraken-dir", required=True,
                        help="directory containing per-sample subdirs with *_kraken2.report")
    parser.add_argument("--out-svg",    required=True, help="output SVG path")
    parser.add_argument("--threshold",  type=float, default=1.0,
                        help="contamination threshold %% for highlighting (default: 1.0)")
    args = parser.parse_args()

    df = collect_data(args.kraken_dir)
    print(f"Samples: {len(df)}", file=sys.stderr)
    print(f"Max contamination: {df['human_pct'].max():.2f}% ({df.loc[df['human_pct'].idxmax(), 'sample']})",
          file=sys.stderr)
    print(f"Above threshold ({args.threshold}%): "
          f"{(df['human_pct'] >= args.threshold).sum()} samples", file=sys.stderr)

    plot(df, args.out_svg, args.threshold)

    # Also write TSV summary
    tsv_path = args.out_svg.replace(".svg", ".tsv")
    df.sort_values("human_pct", ascending=False).to_csv(tsv_path, sep="\t", index=False)
    print(f"Written → {tsv_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
