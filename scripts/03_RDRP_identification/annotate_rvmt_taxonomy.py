#!/usr/bin/env python3
"""
Join DIAMOND blastp hits (query vs. RVMT RCR90 RdRp reference set) with RVMT
taxonomy and add a %identity-based confidence tier.

Input diamond outfmt 6 columns:
    qseqid sseqid pident length evalue bitscore qcovhsp scovhsp
RVMT info tsv (RiboV1.6_Info.tsv) columns used (0-indexed):
    col2  (RCR90)                     = sseqid join key
    col15-19 (Phylum..Genus)          = taxonomy

Output TSV columns:
    qseqid sseqid hit_rank pident length evalue bitscore qcovhsp scovhsp
    Phylum Class Order Family Genus confidence_tier

hit_rank is 1..N per query, ranked by bitscore descending (ties broken by
evalue ascending), preserving every hit DIAMOND reported (no filtering beyond
what --evalue/--max-target-seqs already applied at search time) so a curator
can fall back to hit_rank 2+ when the top hit's taxonomy is uninformative.

Usage:
    python annotate_rvmt_taxonomy.py \
        --diamond  hits.tsv \
        --info-tsv RiboV1.6_Info.tsv \
        --output   hits.annotated.tsv
"""

import argparse
import sys

import pandas as pd

DIAMOND_COLS = ["qseqid", "sseqid", "pident", "length", "evalue",
                 "bitscore", "qcovhsp", "scovhsp"]


TAXA_RANKS = ["Phylum", "Class", "Order", "Family", "Genus"]


def load_taxonomy(info_tsv_path):
    """RCR90 tip label (col3, 1-indexed) -> [Phylum, Class, Order, Family, Genus].

    An RCR90 ID can appear more than once (e.g. once as a "Lvl 0 - Megatree
    leaves." row and again as a "Lvl 1 - BLASTp match ..." row referencing the
    same tip); prefer the Lvl 0 row since it's the curated tree-tip entry.
    """
    lineage = {}
    lvl0 = set()
    with open(info_tsv_path) as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 20:
                continue
            tip = parts[2].strip()
            if not tip:
                continue
            is_lvl0 = parts[11].strip().startswith("Lvl 0")
            if tip in lineage and tip in lvl0 and not is_lvl0:
                continue
            lineage[tip] = [p.strip() for p in parts[15:20]]
            if is_lvl0:
                lvl0.add(tip)
    return lineage


def confidence_tier(pident):
    if pident >= 50:
        return "high (species/genus-level)"
    if pident >= 30:
        return "medium (family/order-level)"
    return "low (weak homology)"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diamond",  required=True)
    parser.add_argument("--info-tsv", required=True)
    parser.add_argument("--output",   required=True)
    args = parser.parse_args()

    print("Loading diamond hits...", file=sys.stderr)
    df = pd.read_csv(args.diamond, sep="\t", header=None, names=DIAMOND_COLS, dtype=str)
    if df.empty:
        pd.DataFrame(columns=DIAMOND_COLS + ["hit_rank"] + TAXA_RANKS + ["confidence_tier"]) \
            .to_csv(args.output, sep="\t", index=False)
        print("No hits — empty output.", file=sys.stderr)
        return

    for col in ["pident", "length", "evalue", "bitscore", "qcovhsp", "scovhsp"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    print(f"  {len(df)} hits for {df['qseqid'].nunique()} queries.", file=sys.stderr)

    df = df.sort_values(["qseqid", "bitscore", "evalue"], ascending=[True, False, True])
    df["hit_rank"] = df.groupby("qseqid").cumcount() + 1

    print("Loading RVMT taxonomy...", file=sys.stderr)
    lineage_map = load_taxonomy(args.info_tsv)
    print(f"  {len(lineage_map)} RVMT sequences with taxonomy.", file=sys.stderr)

    taxa = df["sseqid"].map(lineage_map)
    for i, rank in enumerate(TAXA_RANKS):
        df[rank] = taxa.apply(lambda t, i=i: t[i] if isinstance(t, list) else "")
    df["confidence_tier"] = df["pident"].apply(confidence_tier)

    out_cols = DIAMOND_COLS + ["hit_rank"] + TAXA_RANKS + ["confidence_tier"]
    df = df[out_cols]
    df.to_csv(args.output, sep="\t", index=False)

    n_no_tax = (df["Phylum"] == "").sum()
    print(f"Wrote {len(df)} annotated hits -> {args.output}", file=sys.stderr)
    if n_no_tax:
        print(f"  Warning: {n_no_tax} hits had no RVMT taxonomy match.", file=sys.stderr)


if __name__ == "__main__":
    main()
