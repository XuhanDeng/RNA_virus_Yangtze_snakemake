#!/usr/bin/env python3
"""
Generate iTOL DATASET_COLORSTRIP files for ICTV reference tips colored by taxonomy.

One file per rank: Order, Family, Genus.
Only ICTV tips (prefix ICTV_) get a color; contigs and ESvirtu tips are skipped.

Inputs:
  --vmr         VMR_latest.xlsx
  --vmr-sheet   sheet name (e.g. "VMR MSL41")
  --lookup      my_tree.renamed_lookup.tsv
  --outdir      output directory
"""

import argparse
import os
import re
import sys

try:
    import openpyxl
except ImportError:
    sys.exit("ERROR: openpyxl not installed.")

RANKS = ["Phylum", "Class", "Order", "Family", "Genus"]

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


def assign_colors(values):
    """Assign a color to each distinct value (sorted for reproducibility)."""
    return {v: _PALETTE[i % len(_PALETTE)] for i, v in enumerate(sorted(values))}


def load_vmr(vmr_path, vmr_sheet):
    wb = openpyxl.load_workbook(vmr_path, read_only=True)
    if vmr_sheet and vmr_sheet in wb.sheetnames:
        ws = wb[vmr_sheet]
    else:
        ws = next((wb[s] for s in wb.sheetnames if "VMR" in s or "MSL" in s), wb.active)

    rows = list(ws.iter_rows(values_only=True))
    header = list(rows[0])
    acc_col = header.index("Virus GENBANK accession")
    rank_cols = {rank: header.index(rank) for rank in RANKS}

    acc2tax = {}
    for r in rows[1:]:
        tax = {rank: r[rank_cols[rank]] or "" for rank in RANKS}
        for acc in re.findall(r"[A-Z]+\d+", str(r[acc_col] or "")):
            acc2tax[acc.split(".")[0]] = tax

    print(f"VMR: {len(acc2tax)} accessions loaded", file=sys.stderr)
    return acc2tax


def load_ictv_tips(lookup_path):
    tips = []
    with open(lookup_path) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2 and parts[1].startswith("ICTV_"):
                tips.append(parts[1])
    print(f"ICTV tips: {len(tips)}", file=sys.stderr)
    return tips


def write_colorstrip(outpath, rank, ictv_tips, acc2tax):
    # Collect values that actually appear
    tip_values = {}
    for tip in ictv_tips:
        acc = tip[5:].split(".")[0]
        val = acc2tax.get(acc, {}).get(rank, "")
        tip_values[tip] = val

    used_values = sorted({v for v in tip_values.values() if v})
    colors = assign_colors(used_values)

    with open(outpath, "w") as out:
        out.write("DATASET_COLORSTRIP\n")
        out.write("SEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\tICTV {rank}\n")
        out.write("COLOR\t#000000\n")
        out.write(f"STRIP_WIDTH\t50\n")
        out.write(f"MARGIN\t5\n")
        out.write(f"LEGEND_TITLE\t{rank}\n")
        out.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(used_values)) + "\n")
        out.write("LEGEND_COLORS\t" + "\t".join(colors[v] for v in used_values) + "\n")
        out.write("LEGEND_LABELS\t" + "\t".join(used_values) + "\n")
        out.write("DATA\n")
        for tip, val in tip_values.items():
            if val:
                out.write(f"{tip}\t{colors[val]}\t{val}\n")

    print(f"Written {sum(1 for v in tip_values.values() if v)} entries "
          f"({len(used_values)} {rank}s) → {outpath}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vmr",        required=True)
    parser.add_argument("--vmr-sheet",  default=None)
    parser.add_argument("--lookup",     required=True)
    parser.add_argument("--outdir",     required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    acc2tax   = load_vmr(args.vmr, getattr(args, "vmr_sheet", None))
    ictv_tips = load_ictv_tips(args.lookup)

    for rank in RANKS:
        outpath = os.path.join(args.outdir, f"itol_ictv_{rank.lower()}.txt")
        write_colorstrip(outpath, rank, ictv_tips, acc2tax)


if __name__ == "__main__":
    main()
