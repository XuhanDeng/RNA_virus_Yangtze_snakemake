#!/usr/bin/env python3
"""
Bar chart of phylum-level assignment from gappa per_query phylum_annotation.tsv.
X-axis: phylum; Y-axis: number of queries.
Bars are coloured consistently with the iTOL colorstrip palette.
Unassigned queries (if any) are shown in grey at the right.

Usage:
    python plot_gappa_phylum.py --input phylum_annotation.tsv --out figure.svg [--title "..."]
"""

import argparse
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


_PALETTE = [
    "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
    "#42D4F4", "#F032E6", "#BFEF45", "#FABED4", "#469990",
    "#DCBEFF", "#9A6324", "#FFFAC8", "#800000", "#AAFFC3",
    "#808000", "#FFD8B1", "#000075", "#A9A9A9", "#FF4500",
    "#00CED1", "#FF1493", "#7FFF00", "#DC143C", "#00BFFF",
]
_UNASSIGNED_COLOR = "#CCCCCC"


def load_phylum(path: str) -> dict:
    counts = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            p = row["phylum"].strip() or "Unassigned"
            counts[p] = counts.get(p, 0) + 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",  required=True, help="phylum_annotation.tsv")
    parser.add_argument("--out",    required=True, help="output SVG path")
    parser.add_argument("--title",  default="Phylum Assignment (gappa EPA-ng)",
                        help="figure title")
    args = parser.parse_args()

    counts = load_phylum(args.input)
    total  = sum(counts.values())

    # Sort by count descending; Unassigned always last
    named = sorted(((v, k) for k, v in counts.items() if k != "Unassigned"), reverse=True)
    phyla  = [k for _, k in named]
    values = [v for v, _ in named]
    if "Unassigned" in counts:
        phyla.append("Unassigned")
        values.append(counts["Unassigned"])

    # Assign colours (same order as iTOL script: sorted alphabetically)
    named_sorted = sorted(k for k in counts if k != "Unassigned")
    palette_map  = {p: _PALETTE[i % len(_PALETTE)] for i, p in enumerate(named_sorted)}
    palette_map["Unassigned"] = _UNASSIGNED_COLOR
    colors = [palette_map[p] for p in phyla]

    fig, ax = plt.subplots(figsize=(max(5, len(phyla) * 1.2), 5))
    bars = ax.bar(phyla, values, color=colors, width=0.6, edgecolor="white", linewidth=0.5)

    # Value labels
    for bar, val in zip(bars, values):
        pct = val / total * 100
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + total * 0.005,
                f"{val}\n({pct:.1f}%)",
                ha="center", va="bottom", fontsize=8)

    ax.set_ylabel("Number of queries", fontsize=10)
    ax.set_title(f"{args.title}\n(n = {total:,})", fontsize=11, fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    plt.xticks(rotation=30, ha="right", fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, max(values) * 1.18)

    plt.tight_layout()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Written → {args.out}", file=sys.stderr)

    tsv = args.out.replace(".svg", ".tsv")
    with open(tsv, "w") as fh:
        fh.write("phylum\tcount\tpercent\n")
        for p, v in zip(phyla, values):
            fh.write(f"{p}\t{v}\t{v/total*100:.2f}\n")
    print(f"Written → {tsv}", file=sys.stderr)
    plt.close()


if __name__ == "__main__":
    main()
