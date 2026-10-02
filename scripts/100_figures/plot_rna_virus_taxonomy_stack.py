#!/usr/bin/env python3
"""
Stacked bar plot of per-sample RNA-virus TPM abundance, grouped by a single
unified taxonomic level (tax_phylum / tax_class / tax_order / tax_family /
tax_genus) from rna_virus_recalc_tpm.annotated.tsv.

Rows with no taxonomy call at the given level (tax_<level> is NA -- e.g.
unknown contigs with no RdRp/RVMT hit) are grouped into a single "Unknown"
segment, so each bar still sums to that sample's full RNA-virus TPM total.

To keep the legend readable, only the top-N taxa by total TPM (across all
samples) are shown individually; everything else is collapsed into "Other".
"Unknown" is never collapsed into "Other" -- it is always its own segment.
"""

import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import pandas as pd

UNKNOWN_LABEL = "Unknown"
OTHER_LABEL = "Other"


def sample_sort_key(count_col: str):
    """Sort by the trailing numeric sample index (e.g. ..._42_tpm -> 42)."""
    match = re.search(r"(\d+)_tpm$", count_col)
    return int(match.group(1)) if match else float("inf")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="rna_virus_recalc_tpm.annotated.tsv")
    parser.add_argument("--tax-column", required=True, help="e.g. tax_phylum")
    parser.add_argument("--top-n", type=int, default=12,
                         help="max distinct taxa shown individually; rest collapsed into 'Other'")
    parser.add_argument("--output", required=True, help="output figure path (.svg/.png)")
    parser.add_argument("--output-table", required=True,
                         help="output CSV path with per-sample, per-taxon TPM and relative abundance")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t", low_memory=False)
    if args.tax_column not in df.columns:
        sys.exit(f"ERROR: column '{args.tax_column}' not found in {args.input}")

    tpm_cols = sorted([c for c in df.columns if c.endswith("_tpm")], key=sample_sort_key)
    if not tpm_cols:
        sys.exit("No *_tpm columns found in input.")

    df["_taxon"] = df[args.tax_column].fillna(UNKNOWN_LABEL)

    # sum TPM per (sample, taxon)
    grouped = df.groupby("_taxon")[tpm_cols].sum(min_count=1).fillna(0)

    # rank taxa (excluding Unknown) by total TPM across all samples
    non_unknown = grouped.drop(index=UNKNOWN_LABEL, errors="ignore")
    ranked = non_unknown.sum(axis=1).sort_values(ascending=False)
    top_taxa = list(ranked.index[:args.top_n])
    other_taxa = list(ranked.index[args.top_n:])

    if other_taxa:
        other_sum = grouped.loc[other_taxa].sum(axis=0)
        grouped = grouped.drop(index=other_taxa)
        grouped.loc[OTHER_LABEL] = other_sum

    # order rows: top taxa by abundance first, then Other, then Unknown last
    order = top_taxa + ([OTHER_LABEL] if other_taxa else []) + \
            ([UNKNOWN_LABEL] if UNKNOWN_LABEL in grouped.index else [])
    grouped = grouped.reindex(order, fill_value=0)

    # percent of each sample's combined RNA-virus TPM total
    totals = grouped.sum(axis=0)
    totals_safe = totals.replace(0, 1)
    pct = grouped.div(totals_safe, axis=1) * 100

    samples = [c.removesuffix("_tpm") for c in tpm_cols]
    grouped_named = grouped.copy()
    grouped_named.columns = samples
    pct.columns = samples

    # ── CSV table: long format, one row per (sample, taxon) ────────────────────
    table = pd.concat(
        [
            grouped_named.stack().rename("tpm"),
            pct.stack().rename("relative_abundance_pct"),
        ],
        axis=1,
    ).reset_index()
    table.columns = [args.tax_column, "sample", "tpm", "relative_abundance_pct"]
    table = table[["sample", args.tax_column, "tpm", "relative_abundance_pct"]]
    table.to_csv(args.output_table, index=False)
    print(f"Wrote {args.output_table}", file=sys.stderr)

    # ── plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(10, len(samples) * 0.35), 7))

    cmap = plt.get_cmap("tab20")
    colors = {taxon: cmap(i % 20) for i, taxon in enumerate(pct.index) if taxon not in (OTHER_LABEL, UNKNOWN_LABEL)}
    colors[OTHER_LABEL] = "#999999"
    colors[UNKNOWN_LABEL] = "#CCCCCC"

    bottom = pd.Series(0.0, index=samples)
    for taxon in pct.index:
        values = pct.loc[taxon]
        ax.bar(samples, values, bottom=bottom, label=taxon, color=colors[taxon])
        bottom = bottom + values

    ax.set_ylabel("Percent of RNA-virus TPM (%)")
    ax.set_xlabel("Sample")
    ax.set_title(f"RNA-virus TPM composition per sample ({args.tax_column})")
    ax.set_ylim(0, 100)
    ax.tick_params(axis="x", rotation=90)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, title=args.tax_column)

    fig.tight_layout()
    fig.savefig(args.output)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
