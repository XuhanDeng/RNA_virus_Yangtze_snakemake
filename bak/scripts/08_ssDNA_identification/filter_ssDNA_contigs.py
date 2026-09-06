#!/usr/bin/env python3
"""
Filter GeNomad virus_summary.tsv to ssDNA virus candidates.

Criteria (all must pass):
  1. taxonomy contains at least one keyword from --ssDNA-keywords
  2. topology != --exclude-topology  (default: Provirus)
  3. n_hallmarks >= --min-hallmarks  (default: 1)

Outputs:
  --out-list   : plain text list of passing seq_names (one per line, no header)
  --out-tsv    : filtered summary rows (TSV with header)
"""

import argparse
import sys
import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--summary",           required=True,
                   help="GeNomad *_virus_summary.tsv")
    p.add_argument("--out-list",          required=True,
                   help="Output: one seq_name per line")
    p.add_argument("--out-tsv",           required=True,
                   help="Output: filtered TSV rows")
    p.add_argument("--ssDNA-keywords",    nargs="+", required=True,
                   help="Taxonomy substrings that mark ssDNA viruses")
    p.add_argument("--exclude-topology",  default="Provirus",
                   help="Topology value to exclude (default: Provirus)")
    p.add_argument("--min-hallmarks",     type=int, default=1,
                   help="Minimum n_hallmarks (default: 1)")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    df = pd.read_csv(args.summary, sep="\t", dtype=str)

    if df.empty:
        print(f"[filter_ssDNA_contigs] {args.summary}: empty summary, writing empty outputs",
              file=sys.stderr)
        df.to_csv(args.out_tsv, sep="\t", index=False)
        open(args.out_list, "w").close()
        return

    # coerce numeric columns
    df["n_hallmarks"] = pd.to_numeric(df["n_hallmarks"], errors="coerce").fillna(0)

    # filter topology
    mask = df["topology"] != args.exclude_topology

    # filter hallmarks
    mask &= df["n_hallmarks"] >= args.min_hallmarks

    # filter taxonomy — keep row if ANY keyword appears in the taxonomy string
    tax_col = df["taxonomy"].fillna("")
    keywords = args.ssDNA_keywords
    tax_mask = tax_col.apply(lambda t: any(kw in t for kw in keywords))
    mask &= tax_mask

    filtered = df[mask].copy()

    print(
        f"[filter_ssDNA_contigs] {args.summary}: "
        f"{len(df)} total → {len(filtered)} ssDNA candidates "
        f"(topology!={args.exclude_topology}, hallmarks>={args.min_hallmarks}, "
        f"keywords={keywords})",
        file=sys.stderr,
    )

    filtered.to_csv(args.out_tsv, sep="\t", index=False)
    filtered["seq_name"].to_csv(args.out_list, index=False, header=False)


if __name__ == "__main__":
    main()
