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
  bait_category      — tier (tier1a/tier1b/tier2) or "RVMT"
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


def load_tier_summary(path: str) -> dict[str, str]:
    """Return {contig_id: tier} for tier1a/1b/2 entries in tier_summary.tsv."""
    contig_to_tier = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            tier_col   = header.index("tier")
            contig_col = header.index("contig")
        except ValueError as e:
            sys.exit(f"tier_summary.tsv missing expected column: {e}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(tier_col, contig_col):
                continue
            tier   = parts[tier_col]
            contig = parts[contig_col]
            if tier in ("tier1a", "tier1b", "tier2") and contig:
                contig_to_tier.setdefault(contig, tier)
    print(f"RdRp bait contigs loaded from tier_summary: {len(contig_to_tier)}", file=sys.stderr)
    return contig_to_tier


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",        required=True)
    parser.add_argument("--output",       required=True,
                        help="All passing hits TSV")
    parser.add_argument("--bait-map",     required=True,
                        help="Per-target bait mapping summary TSV")
    parser.add_argument("--tier-summary", required=True,
                        help="tier_summary.tsv — used to label bait source and category")
    parser.add_argument("--evalue",       type=float, default=1e-9)
    parser.add_argument("--min-pident",   type=float, default=95.0)
    parser.add_argument("--min-tcov",     type=float, default=0.95)
    args = parser.parse_args()

    contig_to_tier = load_tier_summary(args.tier_summary)

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
        lambda q: "own_rdrp" if q in contig_to_tier else "RVMT"
    )
    filtered["bait_category"] = filtered["query"].apply(
        lambda q: contig_to_tier.get(q, "RVMT")
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


if __name__ == "__main__":
    main()
