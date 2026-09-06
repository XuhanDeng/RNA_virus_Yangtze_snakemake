#!/usr/bin/env python3
"""
Filter mmseqs2 convertalis output using stringent thresholds, and produce a
bait-to-target mapping table.

Columns expected (from --format-output):
  query, target, pident, alnlen, mismatch, gapopen,
  qstart, qend, tstart, tend, evalue, bits, qlen, tlen, tcov

Filters applied (from Zayed et al. 2022):
  E-value   < 1e-9
  pident    > 95%
  tcov      >= 95%  (target coverage — envelopment criterion)

bait_target_map.tsv columns:
  recovered_contig   — target contig ID
  bait_source        — "own_rdrp" or "RVMT"
  bait_category      — "RdRpCATCH"/"LucaProt" (from all_samples_rdrp_merged.tsv) or "RVMT"
  best_bait_contig   — query ID of the lowest-evalue bait hit
  pident             — percent identity of that hit
  tcov               — target coverage of that hit
  evalue             — e-value of that hit
  n_bait_hits        — total distinct bait contigs that hit this target
"""

import argparse
import sys

import pandas as pd


COLS = [
    "query", "target", "pident", "alnlen", "mismatch", "gapopen",
    "qstart", "qend", "tstart", "tend", "evalue", "bits", "qlen", "tlen", "tcov",
]


def load_rdrp_merged(path: str) -> dict[str, str]:
    """Return {contig_id: source} ("RdRpCATCH"/"LucaProt") from all_samples_rdrp_merged.tsv."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    try:
        contig_to_source = dict(zip(df["contig"], df["source"]))
    except KeyError as e:
        sys.exit(f"all_samples_rdrp_merged.tsv missing expected column: {e}")
    print(f"RdRp bait contigs loaded from rdrp_merged: {len(contig_to_source)}", file=sys.stderr)
    return contig_to_source


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",        required=True)
    parser.add_argument("--output",       required=True,
                        help="All passing hits TSV")
    parser.add_argument("--bait-map",     required=True,
                        help="Per-target bait mapping summary TSV")
    parser.add_argument("--bait-venn",    required=True,
                        help="Venn summary TSV: own_rdrp only / RVMT only / both")
    parser.add_argument("--rdrp-merged",  required=True,
                        help="all_samples_rdrp_merged.tsv — used to label bait source and category")
    parser.add_argument("--evalue",       type=float, default=1e-9)
    parser.add_argument("--min-pident",   type=float, default=95.0)
    parser.add_argument("--min-tcov",     type=float, default=0.95)
    args = parser.parse_args()

    contig_to_source = load_rdrp_merged(args.rdrp_merged)

    df = pd.read_csv(args.input, sep="\t", header=None, names=COLS, dtype=str)
    print(f"Total hits before filter: {len(df)}", file=sys.stderr)

    df["pident"] = pd.to_numeric(df["pident"], errors="coerce")
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["tcov"]   = pd.to_numeric(df["tcov"],   errors="coerce")

    mask = (
        (df["evalue"]  <  args.evalue)     &
        (df["pident"]  >  args.min_pident) &
        (df["tcov"]    >= args.min_tcov)
    )
    filtered = df[mask].copy()
    print(
        f"After filter (evalue<{args.evalue}, pident>{args.min_pident}, "
        f"tcov>={args.min_tcov}): {len(filtered)} hits, "
        f"{filtered['target'].nunique()} unique recovered contigs",
        file=sys.stderr,
    )

    filtered.to_csv(args.output, sep="\t", index=False)

    # Annotate bait source and category
    filtered["bait_source"] = filtered["query"].apply(
        lambda q: "own_rdrp" if q in contig_to_source else "RVMT"
    )
    filtered["bait_category"] = filtered["query"].apply(
        lambda q: contig_to_source.get(q, "RVMT")
    )

    # Best hit = lowest evalue per (target, bait_source)
    best = (
        filtered
        .sort_values("evalue")
        .groupby(["target", "bait_source"], sort=False)
        .first()
        .reset_index()
    )

    # Count distinct bait contigs per target (across all sources)
    n_hits = (
        filtered.groupby("target")["query"]
        .nunique()
        .rename("n_bait_hits")
        .reset_index()
    )

    bait_map = best.merge(n_hits, on="target")
    bait_map = bait_map.rename(columns={
        "target": "recovered_contig",
        "query":  "best_bait_contig",
    })
    bait_map = bait_map[[
        "recovered_contig", "bait_source", "bait_category",
        "best_bait_contig", "pident", "tcov", "evalue", "n_bait_hits",
    ]].sort_values(["bait_source", "evalue"])

    bait_map.to_csv(args.bait_map, sep="\t", index=False)
    print(f"bait_target_map: {len(bait_map)} rows -> {args.bait_map}", file=sys.stderr)

    # Venn summary: which contigs are hit by own_rdrp / RVMT / both
    own_set  = set(filtered[filtered["bait_source"] == "own_rdrp"]["target"])
    rvmt_set = set(filtered[filtered["bait_source"] == "RVMT"]["target"])
    both_set = own_set & rvmt_set
    own_only  = own_set  - both_set
    rvmt_only = rvmt_set - both_set

    venn = pd.DataFrame([
        {"group": "own_rdrp_only", "n_contigs": len(own_only),  "contigs": ";".join(sorted(own_only))},
        {"group": "RVMT_only",     "n_contigs": len(rvmt_only), "contigs": ";".join(sorted(rvmt_only))},
        {"group": "both",          "n_contigs": len(both_set),  "contigs": ";".join(sorted(both_set))},
    ])
    venn.to_csv(args.bait_venn, sep="\t", index=False)
    print(
        f"Venn summary: own_rdrp_only={len(own_only)}, RVMT_only={len(rvmt_only)}, both={len(both_set)}"
        f" -> {args.bait_venn}", file=sys.stderr
    )


if __name__ == "__main__":
    main()
