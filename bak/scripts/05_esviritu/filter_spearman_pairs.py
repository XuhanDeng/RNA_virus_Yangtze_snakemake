#!/usr/bin/env python3
"""
Split a Spearman filtered TSV (seed / target / r / p / r_group) into
known-known and known-unknown pair files, one file per r_group threshold.

Known  : label does NOT start with 'unknown::' or 'acc:'
Unknown: label starts with 'unknown::' or 'acc:'

For known-unknown files an extra column 'contig_name' is added, containing
the clean contig ID with any leading 'unknown::acc:' / 'unknown::' / 'acc:'
prefix stripped.  Example: unknown::acc:AnQing_D_33_0000002790 → AnQing_D_33_0000002790

Outputs (written to --outdir):
  known_known_pair/r{threshold}.tsv
  known_unknown_pair/r{threshold}.tsv
"""

import argparse
import os
import re
import sys

import pandas as pd


def _is_unknown(label: str) -> bool:
    v = str(label).strip()
    vl = v.lower()
    return (
        vl.startswith("unknown::")
        or v.startswith("acc:")
        or vl in ("", "nan", "none", "na", "unknown", "unknow")
    )


def _clean_contig(label: str) -> str:
    """Strip unknown::acc: / unknown:: / acc: prefix to get the bare contig name."""
    v = str(label).strip()
    v = re.sub(r'^unknown::acc:', '', v)
    v = re.sub(r'^unknown::', '', v)
    v = re.sub(r'^acc:', '', v)
    return v


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

    # known-unknown: one side is known, the other is unknown
    df_ku = df[(~seed_unk & target_unk) | (seed_unk & ~target_unk)].copy()

    # add contig_name: clean label of whichever side is unknown
    def _get_contig_name(row):
        if _is_unknown(row["seed"]):
            return _clean_contig(row["seed"])
        return _clean_contig(row["target"])

    df_ku["contig_name"] = df_ku.apply(_get_contig_name, axis=1)

    # ensure known is always seed, unknown is always target
    needs_swap = df_ku["seed"].apply(_is_unknown)
    df_ku.loc[needs_swap, ["seed", "target"]] = \
        df_ku.loc[needs_swap, ["target", "seed"]].values
    df_ku.loc[needs_swap, ["seed_nonzero_samples", "target_nonzero_samples"]] = \
        df_ku.loc[needs_swap, ["target_nonzero_samples", "seed_nonzero_samples"]].values

    # strip unknown::acc: / unknown:: / acc: prefix from seed and target after split
    for col in ("seed", "target"):
        df_ku[col] = df_ku[col].apply(_clean_contig)

    thresholds = sorted(df["r_group"].dropna().unique(), key=float) \
        if "r_group" in df.columns else []

    def _dedup(df_in):
        """Remove duplicate seed-target pairs (A,B) and (B,A) — keep first occurrence."""
        key = df_in.apply(lambda row: tuple(sorted([row["seed"], row["target"]])), axis=1)
        return df_in[~key.duplicated()].reset_index(drop=True)

    for t in thresholds:
        t_str = str(t)
        subset_kk = _dedup(df_kk[df_kk["r_group"] == t])
        subset_ku = _dedup(df_ku[df_ku["r_group"] == t])
        subset_kk.to_csv(os.path.join(kk_dir, f"r{t_str}.tsv"), sep="\t", index=False)
        subset_ku.to_csv(os.path.join(ku_dir, f"r{t_str}.tsv"), sep="\t", index=False)
        print(f"  r>={t_str}: kk={len(subset_kk)}, ku={len(subset_ku)}", file=sys.stderr)

    # ensure all fixed threshold files exist (even if empty)
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
