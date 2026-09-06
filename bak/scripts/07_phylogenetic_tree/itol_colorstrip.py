#!/usr/bin/env python3
"""
Generate two iTOL DATASET_COLORSTRIP files for the FastTree de-novo tree:

  --output       sequence source (rt / esvirtu / assembled)
  --output-tier  RdRp tier (tier1a / tier1b / tier2 / n/a)

Source classification by ID pattern:
  rt        — ID starts with "rt."                              → grey   #808080
  esvirtu   — {Accession}_frame=N or ORF{N}_{Accession}        → orange #FF8C00
  assembled — everything else                                   → blue   #1F78B4

Tier classification from tier_summary.tsv (col1=RdRp_id, col2=tier):
  tier1a    → dark blue  #08519C
  tier1b    → steel blue #2171B5
  tier2     → light blue #6BAED6
  n/a       → light grey #CCCCCC  (RT outgroup, ESvirtu, or not in table)
"""

import argparse
import sys


# ── source colors ────────────────────────────────────────────────────────────
SRC_COLOR = {"rt": "#808080", "esvirtu": "#FF8C00", "assembled": "#1F78B4"}
SRC_LABEL = {"rt": "RT outgroup", "esvirtu": "ESvirtu reference", "assembled": "Assembled contig"}

# ── tier colors ───────────────────────────────────────────────────────────────
TIER_COLOR = {
    "tier1a":    "#D62728",  # red
    "tier1b":    "#FF7F0E",  # orange
    "tier2":     "#2CA02C",  # green
    "esvirtu":   "#9467BD",  # purple
    "rt":        "#7F7F7F",  # grey
    "assembled": "#CCCCCC",  # light grey fallback (assembled but not in tier table)
}
TIER_LABEL = {
    "tier1a":    "Tier 1a (ABCD)",
    "tier1b":    "Tier 1b (depermuted)",
    "tier2":     "Tier 2 (3-motif)",
    "esvirtu":   "ESvirtu reference",
    "rt":        "RT outgroup",
    "assembled": "Unclassified assembled",
}

def classify_source(seq_id: str) -> str:
    if seq_id.startswith("rt."):
        return "rt"
    # ESvirtu sequences are prefixed with "esvirtu_" by workflow 04_03 Step 3
    if seq_id.startswith("esvirtu_"):
        return "esvirtu"
    return "assembled"


def load_tier_summary(path: str) -> dict:
    """Return {RdRp_id: tier} from tier_summary.tsv."""
    id_to_tier = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            id_col   = header.index("RdRp_id")
            tier_col = header.index("tier")
        except ValueError as e:
            sys.exit(f"tier_summary.tsv missing column: {e}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(id_col, tier_col):
                continue
            rdrp_id = parts[id_col].strip().split(":")[0]
            tier    = parts[tier_col].strip()
            if rdrp_id and tier:
                id_to_tier[rdrp_id] = tier
    print(f"Tier entries loaded: {len(id_to_tier)}", file=sys.stderr)
    return id_to_tier


def read_fasta_ids(path: str) -> list:
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0].split(":")[0].strip())
    return ids


def write_colorstrip(path: str, color_map: dict, label_map: dict,
                     rows: list, dataset_label: str, legend_title: str):
    keys = list(color_map.keys())
    with open(path, "w") as out:
        out.write("DATASET_COLORSTRIP\n")
        out.write("SEPARATOR TAB\n")
        out.write(f"DATASET_LABEL\t{dataset_label}\n")
        out.write("COLOR\t#000000\n")
        out.write(f"LEGEND_TITLE\t{legend_title}\n")
        out.write(f"LEGEND_SHAPES\t" + "\t".join(["1"] * len(keys)) + "\n")
        out.write(f"LEGEND_COLORS\t" + "\t".join(color_map[k] for k in keys) + "\n")
        out.write(f"LEGEND_LABELS\t" + "\t".join(label_map[k] for k in keys) + "\n")
        out.write("DATA\n")
        for seq_id, color, label in rows:
            out.write(f"{seq_id}\t{color}\t{label}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--faa",          required=True,
                        help="query_with_rt.faa — all sequences in the FastTree input")
    parser.add_argument("--tier-summary", required=True,
                        help="tier_summary.tsv from motif pipeline")
    parser.add_argument("--output",       required=True,
                        help="Output iTOL file: sequence source")
    parser.add_argument("--output-tier",  required=True,
                        help="Output iTOL file: RdRp tier")
    args = parser.parse_args()

    id_to_tier = load_tier_summary(args.tier_summary)
    seq_ids    = read_fasta_ids(args.faa)

    src_rows  = []
    tier_rows = []
    src_counts  = {}
    tier_counts = {}

    for seq_id in seq_ids:
        # source
        src = classify_source(seq_id)
        src_counts[src] = src_counts.get(src, 0) + 1
        src_rows.append((seq_id, SRC_COLOR[src], SRC_LABEL[src]))

        # tier: assembled sequences get tier from table; RT/ESvirtu use source label
        if seq_id in id_to_tier:
            tier = id_to_tier[seq_id]
        elif src in ("rt", "esvirtu"):
            tier = src
        else:
            # assembled but not found in tier_summary — likely ID mismatch
            print(f"  WARNING: no tier for assembled ID: {seq_id}", file=sys.stderr)
            tier = "assembled"
        tier_counts[tier] = tier_counts.get(tier, 0) + 1
        tier_rows.append((seq_id, TIER_COLOR[tier], TIER_LABEL[tier]))

    print("Source counts:", file=sys.stderr)
    for g in ("rt", "esvirtu", "assembled"):
        print(f"  {SRC_LABEL[g]}: {src_counts.get(g, 0)}", file=sys.stderr)

    print("Tier counts:", file=sys.stderr)
    for t in ("tier1a", "tier1b", "tier2", "esvirtu", "rt"):
        print(f"  {t}: {tier_counts.get(t, 0)}", file=sys.stderr)

    write_colorstrip(args.output, SRC_COLOR, SRC_LABEL, src_rows,
                     "Sequence source", "Sequence source")
    print(f"Written: {args.output}", file=sys.stderr)

    write_colorstrip(args.output_tier, TIER_COLOR, TIER_LABEL, tier_rows,
                     "RdRp tier", "RdRp tier")
    print(f"Written: {args.output_tier}", file=sys.stderr)


if __name__ == "__main__":
    main()
