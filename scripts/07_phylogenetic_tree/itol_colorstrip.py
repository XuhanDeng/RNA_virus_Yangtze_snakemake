#!/usr/bin/env python3
"""
Generate an iTOL DATASET_COLORSTRIP file for the FastTree de-novo tree.

Inputs:
  --faa          query_with_rt.faa — all sequences fed into FastTree
                 (enumerates every tip, including rt. outgroups that are
                 never renamed and therefore absent from the lookup).
  --lookup       renamed_lookup.tsv  (original_id TAB renamed_id)
                 produced by rename_tree_tips.py.
  --rdrp-merged  final_rdrp_merged.tsv  (contig TAB source TAB ...)
  --output       iTOL DATASET_COLORSTRIP file coloured by broad category
                 (Assembled contig / ESvirtu ref / ICTV ref / RT outgroup)
"""

import argparse
import sys

CAT_COLOR = {
    "contig":  "#1F78B4",
    "esvirtu": "#FF7F00",
    "ictv":    "#E31A1C",
    "rt":      "#808080",
}
CAT_LABEL = {
    "contig":  "Assembled contig",
    "esvirtu": "ESvirtu reference",
    "ictv":    "ICTV reference",
    "rt":      "RT outgroup",
}


def read_faa_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def load_lookup(path):
    mapping = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                mapping[parts[0]] = parts[1]
    return mapping


def load_merged(path):
    contig_to_src = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            c_col = header.index("contig")
            s_col = header.index("source")
        except ValueError as e:
            sys.exit(f"final_rdrp_merged.tsv missing column: {e}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(c_col, s_col):
                contig_to_src[parts[c_col]] = parts[s_col]
    print(f"Merged entries loaded: {len(contig_to_src)}", file=sys.stderr)
    return contig_to_src


def classify_tip(tip, contig_to_src):
    if tip.startswith("rt."):
        return "rt"
    if tip.startswith("esvirtu_"):
        return "esvirtu"
    if tip.startswith("ICTV_"):
        return "ictv"
    if tip not in contig_to_src:
        print(f"  WARNING: no source found for tip: {tip}", file=sys.stderr)
    return "contig"


def write_colorstrip(path, color_map, label_map, rows, dataset_label, legend_title):
    keys = list(color_map.keys())
    with open(path, "w") as out:
        out.write("DATASET_COLORSTRIP\n")
        out.write("SEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\t{dataset_label}\n")
        out.write("COLOR\t#000000\n")
        out.write(f"LEGEND_TITLE\t{legend_title}\n")
        out.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(keys)) + "\n")
        out.write("LEGEND_COLORS\t" + "\t".join(color_map[k] for k in keys) + "\n")
        out.write("LEGEND_LABELS\t" + "\t".join(label_map[k] for k in keys) + "\n")
        out.write("DATA\n")
        for seq_id, color, label in rows:
            out.write(f"{seq_id}\t{color}\t{label}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--faa",          required=True,
                        help="query_with_rt.faa fed into FastTree")
    parser.add_argument("--lookup",       required=True,
                        help="renamed_lookup.tsv from rename_tree_tips.py")
    parser.add_argument("--rdrp-merged",  required=True,
                        help="final_rdrp_merged.tsv (contig + source columns)")
    parser.add_argument("--output",       required=True,
                        help="output iTOL DATASET_COLORSTRIP file")
    args = parser.parse_args()

    faa_ids       = read_faa_ids(args.faa)
    lookup        = load_lookup(args.lookup)
    contig_to_src = load_merged(args.rdrp_merged)

    cat_rows   = []
    cat_counts = {}

    for orig_id in faa_ids:
        # Newick truncates LucaProt IDs at the first ':' (branch-length delimiter)
        # so ORF2_Sample_contig:start:end becomes ORF2_Sample_contig in the tree
        tree_id = orig_id.split(":")[0]
        # use renamed ID if it was renamed, otherwise keep the truncated tree ID
        tip = lookup.get(tree_id, tree_id)
        cat = classify_tip(tip, contig_to_src)
        cat_counts[cat] = cat_counts.get(cat, 0) + 1
        cat_rows.append((tip, CAT_COLOR[cat], CAT_LABEL[cat]))

    print("Category counts:", file=sys.stderr)
    for c in CAT_COLOR:
        print(f"  {CAT_LABEL[c]}: {cat_counts.get(c, 0)}", file=sys.stderr)

    write_colorstrip(args.output, CAT_COLOR, CAT_LABEL, cat_rows,
                     "Sequence category", "Sequence category")
    print(f"Written: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
