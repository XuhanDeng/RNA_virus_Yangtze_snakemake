#!/usr/bin/env python3
"""
Generate iTOL DATASET_CONNECTION files from Spearman correlation tables.

Inputs:
  --esvirtu-table   esvirtu_info_table.tsv  (subspecies + Accession columns)
  --lookup          my_tree.renamed_lookup.tsv  (original_id TAB renamed_id)
  --known-known-dir directory with r0.6.tsv … r0.9.tsv (esvirtu ↔ esvirtu)
  --known-unknown-dir directory with r0.6.tsv … r0.9.tsv (esvirtu ↔ contig)
  --outdir          output directory for 4 iTOL files

One output file per threshold (r0.6 … r0.9), each containing all pairs from
that threshold file.  Lines are coloured by 0.1-bin of the actual ρ value:
  [0.6, 0.7)  →  #4DAF4A  (green)
  [0.7, 0.8)  →  #FF7F00  (orange)
  [0.8, 0.9)  →  #E31A1C  (red)
  [0.9, 1.0]  →  #984EA3  (purple)
"""

import argparse
import os
import sys
from typing import Optional

RANK_COLOR = {
    0.6: "#4DAF4A",
    0.7: "#FF7F00",
    0.8: "#E31A1C",
    0.9: "#984EA3",
}
RANK_LABEL = {
    0.6: "r 0.6-0.7",
    0.7: "r 0.7-0.8",
    0.8: "r 0.8-0.9",
    0.9: "r ≥ 0.9",
}
THRESHOLDS = [0.6, 0.7, 0.8, 0.9]


def r_bin(r: float) -> float:
    """Return the lower bound of the 0.1-bin for r."""
    for t in reversed(THRESHOLDS):
        if r >= t:
            return t
    return THRESHOLDS[0]


def load_lookup(path: str) -> dict:
    """Return {original_id: renamed_id} from rename_tree_tips lookup TSV."""
    mapping = {}
    with open(path) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                mapping[parts[0]] = parts[1]
    return mapping


def build_esvirtu_name_map(esvirtu_table: str, lookup: dict) -> dict:
    """
    Return {t__SpeciesName: renamed_tip} by resolving accessions through lookup.
    For multi-accession entries, use the first accession that has a tree tip.
    """
    # Build accession -> renamed_tip from lookup values starting with esvirtu_
    acc_to_tip = {}
    for orig, renamed in lookup.items():
        if orig.startswith("esvirtu_"):
            # orig = esvirtu_AF055887.1_frame=...
            acc = orig[len("esvirtu_"):].split("_frame=")[0]
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


def build_contig_tip_map(lookup: dict) -> dict:
    """
    Return {contig_id: renamed_tip} for non-esvirtu, non-ICTV, non-rt entries.
    Contig renamed tips are the values from lookup where original contains _frame=
    and does not start with esvirtu_ or ICTV_.
    """
    contig_map = {}
    for orig, renamed in lookup.items():
        if orig.startswith("esvirtu_") or orig.startswith("ICTV_") or orig.startswith("rt."):
            continue
        if "_frame=" not in orig:
            continue
        # contig_id is the renamed tip itself (e.g. DongtingLake_D_23_0000000804)
        contig_map[renamed] = renamed
    return contig_map


def resolve_tip(name: str, esvirtu_name_map: dict, contig_tip_map: dict,
                lookup_by_renamed: dict) -> Optional[str]:
    """Resolve a seed/target name to a tree tip ID, or None if not found."""
    # ESvirtu species name
    if name.startswith("t__"):
        return esvirtu_name_map.get(name)
    # Contig ID — appears directly as a renamed tip
    if name in contig_tip_map:
        return name
    return None


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


def write_connection(outpath: str, rows: list, threshold: float):
    """Write iTOL DATASET_CONNECTION file."""
    colors = sorted(RANK_COLOR.keys())
    active = [t for t in colors if t >= threshold]

    with open(outpath, "w") as out:
        out.write("DATASET_CONNECTION\n")
        out.write("SEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\tCorrelation r≥{threshold}\n")
        out.write("COLOR\t#000000\n")
        out.write("DRAW_ARROWS\t0\n")
        out.write("CENTER_CURVES\t1\n")
        out.write("CURVE_ANGLE\t45\n")
        out.write(f"LEGEND_TITLE\tSpearman ρ rank\n")
        out.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(active)) + "\n")
        out.write("LEGEND_COLORS\t" + "\t".join(RANK_COLOR[t] for t in active) + "\n")
        out.write("LEGEND_LABELS\t" + "\t".join(RANK_LABEL[t] for t in active) + "\n")
        out.write("DATA\n")
        for tip1, tip2, r, color in rows:
            out.write(f"{tip1}\t{tip2}\t1\t{color}\tnormal\t{r:.4f}\n")
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

    lookup          = load_lookup(args.lookup)
    esvirtu_name_map = build_esvirtu_name_map(args.esvirtu_table, lookup)
    contig_tip_map  = build_contig_tip_map(lookup)

    # Build reverse lookup: renamed_tip -> True, for fast membership check
    renamed_tips = set(lookup.values())

    missing_seeds   = set()
    missing_targets = set()

    for threshold in THRESHOLDS:
        # input files are bins (e.g. r0.7 = [0.7,0.8)), so read this bin and all higher ones
        bins_to_read = [t for t in THRESHOLDS if t >= threshold]
        rows = []

        for bin_t in bins_to_read:
            fname   = f"r{bin_t}.tsv"
            kk_path = os.path.join(args.known_known_dir,   fname)
            ku_path = os.path.join(args.known_unknown_dir, fname)

            for path, is_ku in [(kk_path, False), (ku_path, True)]:
                if not os.path.exists(path):
                    continue
                pairs = read_pairs(path)
                if not pairs:
                    continue
                for seed, target, r, sn, tn in pairs:
                    if is_ku:
                        if args.min_seed_nonzero   > 0 and sn < args.min_seed_nonzero:
                            continue
                        if args.min_target_nonzero > 0 and tn < args.min_target_nonzero:
                            continue
                    tip1 = resolve_tip(seed,   esvirtu_name_map, contig_tip_map, lookup)
                    tip2 = resolve_tip(target, esvirtu_name_map, contig_tip_map, lookup)
                    if tip1 is None:
                        missing_seeds.add(seed)
                        continue
                    if tip2 is None:
                        missing_targets.add(target)
                        continue
                    color = RANK_COLOR[r_bin(r)]
                    rows.append((tip1, tip2, r, color))

        outpath = os.path.join(args.outdir, f"itol_connection_r{threshold}.txt")
        write_connection(outpath, rows, threshold)

    if missing_seeds:
        print(f"Unresolved seeds ({len(missing_seeds)}): "
              + ", ".join(sorted(missing_seeds)[:10]), file=sys.stderr)
    if missing_targets:
        print(f"Unresolved targets ({len(missing_targets)}): "
              + ", ".join(sorted(missing_targets)[:10]), file=sys.stderr)


if __name__ == "__main__":
    main()
