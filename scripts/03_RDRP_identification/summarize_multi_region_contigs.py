#!/usr/bin/env python3
"""
Summarize contigs that had more than one surviving candidate RdRp region
after overlap dedup (before the longest-per-contig selection), and show
which region the pipeline actually kept.

Reads:
  --regions  all_samples_rdrp_regions.tsv or ICTV_rdrp_regions.tsv
             (all surviving candidates, pre-longest)
  --merged   all_samples_rdrp_merged.tsv or ICTV_rdrp_merged.tsv
             (final, one row per contig)

Works for both the per-sample RNA path (regions table has a `sample` column)
and the single-file ICTV path (no `sample` column -- treated as one group).

Writes:
  --output   one row per (sample, contig, candidate region), with an
             is_kept flag marking which candidate matches the final call.
"""

import argparse
import sys

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regions", required=True)
    parser.add_argument("--merged",  required=True)
    parser.add_argument("--output",  required=True)
    args = parser.parse_args()

    regions = pd.read_csv(args.regions, sep="\t", dtype=str)
    merged  = pd.read_csv(args.merged,  sep="\t", dtype=str)

    if "sample" not in regions.columns:
        regions.insert(0, "sample", "ICTV")

    regions["aa_length_f"] = pd.to_numeric(regions["aa_length"], errors="coerce")

    counts = regions.groupby(["sample", "contig"]).size()
    multi_contigs = set(counts[counts > 1].index)
    print(f"Contigs with >1 surviving region after dedup: {len(multi_contigs)}",
          file=sys.stderr)

    kept_key = {
        (row["contig"], row["source"], row.get("aa_length"))
        for _, row in merged.iterrows()
    }

    rows = []
    for (sample, contig) in sorted(multi_contigs):
        cands = regions[(regions["sample"] == sample) & (regions["contig"] == contig)]
        for _, r in cands.iterrows():
            is_kept = (r["contig"], r["source"], r.get("aa_length")) in kept_key
            rows.append({
                "sample":     sample,
                "contig":     contig,
                "source":     r["source"],
                "aa_length":  r["aa_length"],
                "strand":     r.get("strand"),
                "pssm_ABC":   r.get("pssm_ABC"),
                "is_kept":    is_kept,
            })

    out = pd.DataFrame(rows)
    out.to_csv(args.output, sep="\t", index=False)
    print(f"Wrote {len(out)} rows ({len(multi_contigs)} contigs) -> {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
