#!/usr/bin/env python3
"""
Split a Spearman filtered TSV into known-known and known-unknown pair files.
Adapted from scripts/05_esviritu/filter_spearman_pairs.py for use with
the recalculated RNA-virus-only RPKMF table (rna_virus_recalc_rpkmf.tsv).

Change from original: known/unknown classification uses t__ prefix:
  Known  : label starts with 't__' (ESvirtu reference)
  Unknown: everything else (novel contig bare IDs)

Outputs (written to --outdir):
  known_known_pair/r{threshold}.tsv
  known_unknown_pair/r{threshold}.tsv
"""

import argparse
import os
import sys

import pandas as pd


def _is_unknown(label: str) -> bool:
    v = str(label).strip()
    return not v.startswith("t__")


def _write_empty(path: str, columns: list):
    pd.DataFrame(columns=columns).to_csv(path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Split Spearman filtered TSV into known-known / known-unknown per threshold."
    )
    parser.add_argument("--input",  required=True, help="Spearman filtered TSV")
    parser.add_argument("--outdir", required=True, help="Base output directory")
    args = parser.parse_args()

    kk_dir = os.path.join(args.outdir, "known_known_pair")
    ku_dir = os.path.join(args.outdir, "known_unknown_pair")
    os.makedirs(kk_dir, exist_ok=True)
    os.makedirs(ku_dir, exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str)

    base_cols = ["seed", "target", "r", "p", "r_group",
                 "seed_nonzero_samples", "target_nonzero_samples",
                 "nonzero_intersection", "nonzero_union"]
    ku_cols   = base_cols + ["contig_name"]
    fixed_thresholds = ["0.6", "0.7", "0.8", "0.9"]

    if df.empty or "seed" not in df.columns:
        print(f"Input is empty or missing 'seed' column: {args.input}", file=sys.stderr)
        for t in fixed_thresholds:
            _write_empty(os.path.join(kk_dir, f"r{t}.tsv"), base_cols)
            _write_empty(os.path.join(ku_dir, f"r{t}.tsv"), ku_cols)
        return

    seed_unk   = df["seed"].apply(_is_unknown)
    target_unk = df["target"].apply(_is_unknown)

    df_kk = df[~seed_unk & ~target_unk].copy()
    df_ku = df[(~seed_unk & target_unk) | (seed_unk & ~target_unk)].copy()

    # contig_name: the bare contig ID (already bare in this table, no prefix to strip)
    def _get_contig_name(row):
        if _is_unknown(row["seed"]):
            return str(row["seed"]).strip()
        return str(row["target"]).strip()

    df_ku["contig_name"] = df_ku.apply(_get_contig_name, axis=1)

    # ensure known is always seed, unknown is always target
    needs_swap = df_ku["seed"].apply(_is_unknown)
    df_ku.loc[needs_swap, ["seed", "target"]] = \
        df_ku.loc[needs_swap, ["target", "seed"]].values
    df_ku.loc[needs_swap, ["seed_nonzero_samples", "target_nonzero_samples"]] = \
        df_ku.loc[needs_swap, ["target_nonzero_samples", "seed_nonzero_samples"]].values

    thresholds = sorted(df["r_group"].dropna().unique(), key=float) \
        if "r_group" in df.columns else []

    def _dedup(df_in):
        key = df_in.apply(lambda row: tuple(sorted([row["seed"], row["target"]])), axis=1)
        return df_in[~key.duplicated()].reset_index(drop=True)

    for t in thresholds:
        t_str = str(t)
        subset_kk = _dedup(df_kk[df_kk["r_group"] == t])
        subset_ku = _dedup(df_ku[df_ku["r_group"] == t])
        subset_kk.to_csv(os.path.join(kk_dir, f"r{t_str}.tsv"), sep="\t", index=False)
        subset_ku.to_csv(os.path.join(ku_dir, f"r{t_str}.tsv"), sep="\t", index=False)
        print(f"  r>={t_str}: kk={len(subset_kk)}, ku={len(subset_ku)}", file=sys.stderr)

    for t in fixed_thresholds:
        kk_path = os.path.join(kk_dir, f"r{t}.tsv")
        ku_path = os.path.join(ku_dir, f"r{t}.tsv")
        if not os.path.exists(kk_path):
            _write_empty(kk_path, base_cols)
        if not os.path.exists(ku_path):
            _write_empty(ku_path, ku_cols)

    print(f"Done. kk total={len(df_kk)}, ku total={len(df_ku)}", file=sys.stderr)


if __name__ == "__main__":
    main()
