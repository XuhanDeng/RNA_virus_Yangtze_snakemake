#!/usr/bin/env python3
"""
Generate an iTOL DATASET_COLORSTRIP annotation for the filtered FastTree from
gappa phylum assignments.

Only contig tips are annotated (ESvirtu tips are skipped).
Phylum "Unassigned" gets a light grey color.

Inputs:
  --phylum-tsv   phylum_annotation.tsv (from parse_gappa_per_query.py)
  --out          output iTOL colorstrip file
"""

import argparse
import csv
import os
import re
import sys


_PALETTE = [
    "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
    "#42D4F4", "#F032E6", "#BFEF45", "#FABED4", "#469990",
    "#DCBEFF", "#9A6324", "#FFFAC8", "#800000", "#AAFFC3",
    "#808000", "#FFD8B1", "#000075", "#A9A9A9", "#FF4500",
    "#00CED1", "#FF1493", "#7FFF00", "#DC143C", "#00BFFF",
    "#FF8C00", "#9400D3", "#00FA9A", "#FF6347", "#4682B4",
    "#DA70D6", "#32CD32", "#FF00FF", "#1E90FF", "#FFD700",
    "#8B0000", "#006400", "#00008B", "#8B4513", "#2E8B57",
    "#B8860B", "#4B0082", "#556B2F", "#8B008B", "#008B8B",
    "#CD853F", "#6B8E23",
]
_UNASSIGNED_COLOR = "#CCCCCC"


def rename_tip(tip: str) -> str:
    """Same logic as rename_tree_tips.py — maps original query ID to tree tip label."""
    # LucaProt: ORF{N}_{rest}:{start}:{end}
    m = re.match(r"^ORF\d+_(.+?)(?::\d+:\d+)?$", tip)
    if m:
        return m.group(1)
    # RdRpCATCH: strip _frame=... suffix
    m = re.match(r"^(.+?)_frame=.*$", tip)
    if m:
        return m.group(1)
    return tip


def query_to_tip(query: str) -> str:
    """Strip motif annotation (everything after first space), then rename."""
    bare = query.split(" ")[0]
    return rename_tip(bare)


def assign_colors(phyla: set) -> dict:
    named = sorted(p for p in phyla if p != "Unassigned")
    colors = {p: _PALETTE[i % len(_PALETTE)] for i, p in enumerate(named)}
    colors["Unassigned"] = _UNASSIGNED_COLOR
    return colors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phylum-tsv", required=True,
                        help="phylum_annotation.tsv from parse_gappa_per_query.py")
    parser.add_argument("--out", required=True, help="output iTOL colorstrip .txt file")
    args = parser.parse_args()

    # Load phylum assignments
    tip2phylum = {}
    with open(args.phylum_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tip = query_to_tip(row["query"])
            tip2phylum[tip] = row["phylum"]

    all_phyla = set(tip2phylum.values())
    colors = assign_colors(all_phyla)

    # Legend entries (sorted, Unassigned last)
    legend_phyla = sorted(p for p in all_phyla if p != "Unassigned")
    if "Unassigned" in all_phyla:
        legend_phyla.append("Unassigned")

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("DATASET_COLORSTRIP\n")
        fh.write("SEPARATOR TAB\n")
        fh.write("DATASET_LABEL\tPhylum (gappa)\n")
        fh.write("COLOR\t#888888\n")
        fh.write("COLOR_BRANCHES\t0\n")
        fh.write("STRIP_WIDTH\t40\n")
        fh.write("SHOW_INTERNAL\t0\n")
        fh.write(f"LEGEND_TITLE\tPhylum\n")
        fh.write(f"LEGEND_SHAPES\t" + "\t".join("1" for _ in legend_phyla) + "\n")
        fh.write(f"LEGEND_COLORS\t" + "\t".join(colors[p] for p in legend_phyla) + "\n")
        fh.write(f"LEGEND_LABELS\t" + "\t".join(legend_phyla) + "\n")
        fh.write("DATA\n")
        for tip, phylum in sorted(tip2phylum.items()):
            fh.write(f"{tip}\t{colors[phylum]}\t{phylum}\n")

    print(f"Written → {args.out}  ({len(tip2phylum)} tips, {len(all_phyla)} phyla)",
          file=sys.stderr)
    for p in legend_phyla:
        n = sum(1 for v in tip2phylum.values() if v == p)
        print(f"  {colors[p]}  {p}: {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
