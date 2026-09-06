#!/usr/bin/env python3
"""
Generate iTOL DATASET_CONNECTION files coloured by known virus identity.

Coloring:
  known-known pairs  → black (#000000), both endpoints are esvirtu tips
  known-unknown pairs → color determined by the esvirtu seed species;
                        each unique esvirtu species gets a distinct color

Line style encodes ρ rank (pairs with ρ < 0.7 are skipped):
  [0.7, 0.8)  → dashed
  [0.8, 0.9)  → normal
  [0.9, 1.0]  → normal

All lines share the same width (1).
"""

import argparse
import os
import sys
from typing import Optional


THRESHOLDS = [0.6, 0.7, 0.8, 0.9]


def line_style(r: float) -> Optional[str]:
    """Return iTOL line style for a given ρ, or None to skip."""
    if r >= 0.9:
        return "normal"
    if r >= 0.8:
        return "normal"
    if r >= 0.7:
        return "dashed"
    return None  # skip ρ < 0.7


def load_lookup(path: str) -> dict:
    mapping = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                mapping[parts[0]] = parts[1]
    return mapping


def build_esvirtu_name_map(esvirtu_table: str, lookup: dict) -> dict:
    acc_to_tip = {}
    for orig, renamed in lookup.items():
        # Handles both "esvirtu_ACC_frame=..." and "ORF2_esvirtu_ACC" tip formats
        m = orig.find("esvirtu_")
        if m == -1:
            continue
        tail = orig[m + len("esvirtu_"):]
        acc = tail.split("_frame=")[0]
        acc_to_tip[acc] = renamed

    name_to_tip = {}
    with open(esvirtu_table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        sub_col = header.index("subspecies")
        acc_col = header.index("Accession")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(sub_col, acc_col):
                continue
            name = parts[sub_col]
            accs = parts[acc_col].split(";")
            for acc in accs:
                acc = acc.strip()
                if acc in acc_to_tip:
                    name_to_tip[name] = acc_to_tip[acc]
                    break
    print(f"ESvirtu name→tip mappings: {len(name_to_tip)}", file=sys.stderr)
    return name_to_tip


def build_contig_tip_set(lookup: dict) -> set:
    tips = set()
    for orig, renamed in lookup.items():
        if "esvirtu_" in orig or orig.startswith("ICTV_") or orig.startswith("rt."):
            continue
        if "_frame=" in orig:
            tips.add(renamed)
    return tips


def read_pairs(path: str) -> list:
    """Return list of (seed, target, r, seed_nonzero, target_nonzero)."""
    pairs = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        seed_col   = header.index("seed")
        target_col = header.index("target")
        r_col      = header.index("r")
        sn_col = header.index("seed_nonzero_samples")   if "seed_nonzero_samples"   in header else None
        tn_col = header.index("target_nonzero_samples") if "target_nonzero_samples" in header else None
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(seed_col, target_col, r_col):
                continue
            sn = int(parts[sn_col]) if sn_col is not None and len(parts) > sn_col else 0
            tn = int(parts[tn_col]) if tn_col is not None and len(parts) > tn_col else 0
            pairs.append((parts[seed_col], parts[target_col], float(parts[r_col]), sn, tn))
    return pairs


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


def assign_colors(species_list: list) -> dict:
    """Assign a high-contrast color to each species name."""
    colors = {}
    for i, name in enumerate(sorted(species_list)):
        colors[name] = _PALETTE[i % len(_PALETTE)]
    return colors


def write_connection(outpath: str, rows: list, threshold: float,
                     species_colors: dict, color_to_species: dict):
    # Build legend only from colors that actually appear in rows
    has_black = any(color == "#000000" for _, _, _, color, _ in rows)
    used_species = sorted({
        color_to_species[color]
        for _, _, _, color, _ in rows
        if color != "#000000" and color in color_to_species
    })

    legend_colors = (["#000000"] if has_black else []) + \
                    [species_colors[s] for s in used_species]
    legend_labels = (["Known-Known"] if has_black else []) + \
                    [s.replace("t__", "") for s in used_species]

    with open(outpath, "w") as out:
        out.write("DATASET_CONNECTION\n")
        out.write("SEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\tCorrelation by virus r≥{threshold}\n")
        out.write("COLOR\t#000000\n")
        out.write("DRAW_ARROWS\t0\n")
        out.write("CENTER_CURVES\t1\n")
        out.write("CURVE_ANGLE\t45\n")
        out.write("LEGEND_TITLE\tKnown virus\n")
        if legend_colors:
            out.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_colors)) + "\n")
            out.write("LEGEND_COLORS\t" + "\t".join(legend_colors) + "\n")
            out.write("LEGEND_LABELS\t" + "\t".join(legend_labels) + "\n")
        out.write("DATA\n")
        for tip1, tip2, r, color, style in rows:
            out.write(f"{tip1}\t{tip2}\t1\t{color}\t{style}\t{r:.4f}\n")
    print(f"Written {len(rows)} connections → {outpath}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--esvirtu-table",       required=True)
    parser.add_argument("--lookup",              required=True)
    parser.add_argument("--known-known-dir",     required=True)
    parser.add_argument("--known-unknown-dir",   required=True)
    parser.add_argument("--outdir",              required=True)
    parser.add_argument("--min-seed-nonzero",   type=int, default=0,
                        help="min seed_nonzero_samples for known-unknown pairs (0=no filter)")
    parser.add_argument("--min-target-nonzero", type=int, default=0,
                        help="min target_nonzero_samples for known-unknown pairs (0=no filter)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    lookup           = load_lookup(args.lookup)
    esvirtu_name_map = build_esvirtu_name_map(args.esvirtu_table, lookup)
    contig_tip_set   = build_contig_tip_set(lookup)

    # First pass: collect all esvirtu species that appear as seeds in known-unknown
    seen_species = set()
    for threshold in THRESHOLDS:
        ku_path = os.path.join(args.known_unknown_dir, f"r{threshold}.tsv")
        if not os.path.exists(ku_path):
            continue
        for seed, target, r, sn, tn in read_pairs(ku_path):
            if seed.startswith("t__") and seed in esvirtu_name_map:
                seen_species.add(seed)

    species_colors  = assign_colors(list(seen_species))
    color_to_species = {v: k for k, v in species_colors.items()}
    print(f"Distinct known-virus colors: {len(species_colors)}", file=sys.stderr)

    missing = set()

    for threshold in THRESHOLDS:
        # input files are bins (e.g. r0.7 = [0.7,0.8)), so read this bin and all higher ones
        bins_to_read = [t for t in THRESHOLDS if t >= threshold]
        rows = []

        for bin_t in bins_to_read:
            fname   = f"r{bin_t}.tsv"
            kk_path = os.path.join(args.known_known_dir,   fname)
            ku_path = os.path.join(args.known_unknown_dir, fname)

            # known-known: black (no nonzero filter)
            if os.path.exists(kk_path):
                for seed, target, r, sn, tn in read_pairs(kk_path):
                    style = line_style(r)
                    if style is None:
                        continue
                    tip1 = esvirtu_name_map.get(seed)
                    tip2 = esvirtu_name_map.get(target)
                    if tip1 is None:
                        missing.add(seed)
                        continue
                    if tip2 is None:
                        missing.add(target)
                        continue
                    rows.append((tip1, tip2, r, "#000000", style))

            # known-unknown: color by seed species, apply nonzero filter
            if os.path.exists(ku_path):
                for seed, target, r, sn, tn in read_pairs(ku_path):
                    style = line_style(r)
                    if style is None:
                        continue
                    if args.min_seed_nonzero   > 0 and sn < args.min_seed_nonzero:
                        continue
                    if args.min_target_nonzero > 0 and tn < args.min_target_nonzero:
                        continue
                    tip1 = esvirtu_name_map.get(seed)
                    tip2 = target if target in contig_tip_set else None
                    if tip1 is None:
                        missing.add(seed)
                        continue
                    if tip2 is None:
                        missing.add(target)
                        continue
                    color = species_colors.get(seed, "#888888")
                    rows.append((tip1, tip2, r, color, style))

        outpath = os.path.join(args.outdir, f"itol_connection_r{threshold}.txt")
        write_connection(outpath, rows, threshold, species_colors, color_to_species)

    if missing:
        print(f"Unresolved IDs ({len(missing)}): "
              + ", ".join(sorted(missing)[:10]), file=sys.stderr)


if __name__ == "__main__":
    main()
