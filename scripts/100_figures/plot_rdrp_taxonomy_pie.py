#!/usr/bin/env python3
"""
Pie chart of RdRp cluster-representative taxonomy at a given rank
(rvmt_Phylum / rvmt_Class / rvmt_Order) from
nr_filtered_protein_contig_cluster.tsv, restricted to
is_cluster_representative == True rows.

Usage:
    python plot_rdrp_taxonomy_pie.py --input table.tsv --column rvmt_Phylum \
        --out figure.svg --title "RDRP DATASET"
"""

import argparse
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt


_PALETTE = [
    "#B7E4D8", "#D9D2E9", "#FBE3C2", "#B4D7F0", "#C9E6A3",
    "#F6C7B6", "#E0BBE4", "#FFDFBA", "#BAE1FF", "#D4A5A5",
]
_UNKNOWN_LABEL = "Unknown"
_UNKNOWN_COLOR = "#CCCCCC"


def load_counts(path: str, column: str) -> dict:
    counts = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "is_cluster_representative" not in reader.fieldnames:
            sys.exit(f"ERROR: column 'is_cluster_representative' not found in {path}")
        if column not in reader.fieldnames:
            sys.exit(f"ERROR: column '{column}' not found in {path}")
        for row in reader:
            if row["is_cluster_representative"] != "True":
                continue
            val = (row[column] or "").strip()
            label = val if val else _UNKNOWN_LABEL
            counts[label] = counts.get(label, 0) + 1
    return counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",  required=True, help="nr_filtered_protein_contig_cluster.tsv")
    parser.add_argument("--column", required=True, help="taxonomy column, e.g. rvmt_Phylum")
    parser.add_argument("--out",    required=True, help="output SVG path")
    parser.add_argument("--title",  default="RDRP DATASET", help="figure title")
    args = parser.parse_args()

    counts = load_counts(args.input, args.column)
    total = sum(counts.values())
    if total == 0:
        sys.exit(f"ERROR: no cluster-representative rows found in {args.input}")

    # Sort descending by count; Unknown always last.
    named = sorted(((v, k) for k, v in counts.items() if k != _UNKNOWN_LABEL), reverse=True)
    labels = [k for _, k in named]
    values = [v for v, _ in named]
    if _UNKNOWN_LABEL in counts:
        labels.append(_UNKNOWN_LABEL)
        values.append(counts[_UNKNOWN_LABEL])

    colors = [_PALETTE[i % len(_PALETTE)] for i in range(len(labels))]
    if labels and labels[-1] == _UNKNOWN_LABEL:
        colors[-1] = _UNKNOWN_COLOR

    fig, ax = plt.subplots(figsize=(8, 8))
    wedges, _ = ax.pie(values, colors=colors, startangle=90, counterclock=False,
                        wedgeprops={"edgecolor": "white", "linewidth": 1})

    for wedge, label, val in zip(wedges, labels, values):
        angle = math.radians((wedge.theta2 + wedge.theta1) / 2)
        xy = (wedge.r * math.cos(angle), wedge.r * math.sin(angle))
        xytext = (1.3 * math.cos(angle), 1.3 * math.sin(angle))
        ax.annotate(
            f"{label}\n{val}",
            xy=xy,
            xytext=xytext,
            ha="center", va="center", fontsize=10, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1),
            arrowprops=dict(arrowstyle="-", color="black", lw=0.8),
        )

    ax.set_title(args.title, fontsize=16, fontweight="bold", color="#B8860B", pad=20)

    plt.tight_layout()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    plt.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"Written → {args.out}", file=sys.stderr)

    tsv = args.out.replace(".svg", ".tsv")
    with open(tsv, "w") as fh:
        fh.write(f"{args.column}\tcount\tpercent\n")
        for label, val in zip(labels, values):
            fh.write(f"{label}\t{val}\t{val/total*100:.2f}\n")
    print(f"Written → {tsv}", file=sys.stderr)
    plt.close()


if __name__ == "__main__":
    main()
