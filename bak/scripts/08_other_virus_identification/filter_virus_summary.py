#!/usr/bin/env python3
"""Filter GeNomad virus_summary.tsv: n_hallmarks >= min_hallmarks and topology != exclude_topology."""

import argparse
import sys

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary",           required=True)
    parser.add_argument("--output",            required=True)
    parser.add_argument("--exclude-topology",  default="Provirus")
    parser.add_argument("--min-hallmarks",     type=int, default=1)
    args = parser.parse_args()

    df = pd.read_csv(args.summary, sep="\t", dtype=str)

    if df.empty:
        df.to_csv(args.output, sep="\t", index=False)
        print(f"[filter_virus_summary] {args.summary}: empty, written empty output", file=sys.stderr)
        return

    df["n_hallmarks"] = pd.to_numeric(df["n_hallmarks"], errors="coerce").fillna(0)
    mask = (df["topology"] != args.exclude_topology) & (df["n_hallmarks"] >= args.min_hallmarks)
    filtered = df[mask].copy()

    print(
        f"[filter_virus_summary] {args.summary}: {len(df)} total -> {len(filtered)} passed "
        f"(topology!={args.exclude_topology}, hallmarks>={args.min_hallmarks})",
        file=sys.stderr,
    )

    filtered.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
