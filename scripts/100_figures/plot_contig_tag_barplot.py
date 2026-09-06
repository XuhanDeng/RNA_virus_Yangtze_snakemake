#!/usr/bin/env python3
"""
Stacked bar plot of per-sample read-count percentage by contig_tag.

Combines two sources of RNA-virus read counts:
  - rna_virus_contig_table.tsv: novel contigs tagged by contig_tag
    (tier1a/1b/2/3, bait_tier1a/1b/2/RVMT, genomad_*, unknown).
  - esvirtu_info_table.tsv: ESvirtu reference-mapped hits, tagged by
    genome_type as "esvirtu_<genome_type>" (e.g. esvirtu_ssRNA(+)).

Sums each sample's *_read_count column grouped by tag, converts to percent
of that sample's combined total, and draws one stacked bar per sample.

Rare/mixed genomad_* composite tags (anything other than genomad_ssDNA or
genomad_dsDNA) are collapsed into a single "genomad_other" category to keep
the legend readable.
"""

import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")
# keep SVG text as real <text> elements (editable in Illustrator),
# instead of being converted to paths.
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import pandas as pd

# tags plotted in this fixed order (bottom-to-top in the stack; also legend order).
# ESvirtu categories are stacked at the bottom of each bar.
TAG_ORDER = [
    "esvirtu_ssRNA(+)",
    "identified_RdRp",
    "bait",
    "esvirtu_ssDNA",
    "genomad_ssDNA", "genomad_dsDNA",
    "esvirtu_other", "genomad_other",
    "unknown",
]

KEEP_GENOMAD_TAGS = {"genomad_ssDNA", "genomad_dsDNA"}
KEEP_ESVIRTU_TAGS = {"esvirtu_ssRNA(+)", "esvirtu_ssDNA"}

# RdRpCATCH/LucaProt -> identified_RdRp; all bait_* -> bait
CONTIG_TAG_MERGE = {
    "RdRpCATCH": "identified_RdRp", "LucaProt": "identified_RdRp",
    "bait_RdRpCATCH": "bait", "bait_LucaProt": "bait", "bait_RVMT": "bait",
}


def collapse_genomad_tag(tag: str) -> str:
    if tag.startswith("genomad_") and tag not in KEEP_GENOMAD_TAGS:
        return "genomad_other"
    return tag


def merge_contig_tag(tag: str) -> str:
    return CONTIG_TAG_MERGE.get(tag, tag)


def sample_sort_key(count_col: str):
    """Sort by the trailing numeric sample index (e.g. ..._42_read_count -> 42)."""
    match = re.search(r"(\d+)_read_count$", count_col)
    return int(match.group(1)) if match else float("inf")


def collapse_esvirtu_tag(tag: str) -> str:
    if tag.startswith("esvirtu_") and tag not in KEEP_ESVIRTU_TAGS:
        return "esvirtu_other"
    return tag


def load_contig_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    count_cols = [c for c in df.columns if c.endswith("_read_count")]
    df["tag"] = df["contig_tag"].apply(merge_contig_tag).apply(collapse_genomad_tag)
    return df[["tag"] + count_cols]


def load_esvirtu_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    count_cols = [c for c in df.columns if c.endswith("_read_count")]
    df["tag"] = ("esvirtu_" + df["genome_type"].astype(str)).apply(collapse_esvirtu_tag)
    return df[["tag"] + count_cols]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contig-table",  required=True, help="rna_virus_contig_table.tsv")
    parser.add_argument("--esvirtu-table", required=True, help="esvirtu_info_table.tsv")
    parser.add_argument("--output", required=True, help="output figure path (.svg/.png)")
    parser.add_argument("--output-table", required=True,
                         help="output CSV path with per-sample, per-tag read count and relative abundance")
    parser.add_argument("--include-unknown", action="store_true",
                         help="include contig_tag=='unknown' rows in the percentages")
    args = parser.parse_args()

    contig_df  = load_contig_table(args.contig_table)
    esvirtu_df = load_esvirtu_table(args.esvirtu_table)

    combined = pd.concat([contig_df, esvirtu_df], axis=0, ignore_index=True)

    count_cols = sorted(set(contig_df.columns) & set(esvirtu_df.columns) - {"tag"}, key=sample_sort_key)
    if not count_cols:
        sys.exit("No shared *_read_count columns found between the two input tables.")

    if not args.include_unknown:
        combined = combined[combined["tag"] != "unknown"]

    # sum read counts per (sample, tag)
    grouped = combined.groupby("tag")[count_cols].sum(min_count=1).fillna(0)

    # order rows: fixed TAG_ORDER first, then any leftovers (shouldn't happen)
    present = [t for t in TAG_ORDER if t in grouped.index]
    leftover = [t for t in grouped.index if t not in TAG_ORDER]
    grouped = grouped.reindex(present + leftover, fill_value=0)

    # percent of each sample's combined total
    totals = grouped.sum(axis=0)
    totals_safe = totals.replace(0, 1)
    pct = grouped.div(totals_safe, axis=1) * 100

    samples = [c.removesuffix("_read_count") for c in count_cols]
    grouped_named = grouped.copy()
    grouped_named.columns = samples
    pct.columns = samples

    # ── CSV table: long format, one row per (sample, tag) ───────────────────────
    table = pd.concat(
        [
            grouped_named.stack().rename("read_count"),
            pct.stack().rename("relative_abundance_pct"),
        ],
        axis=1,
    ).reset_index()
    table.columns = ["tag", "sample", "read_count", "relative_abundance_pct"]
    table = table[["sample", "tag", "read_count", "relative_abundance_pct"]]
    table.to_csv(args.output_table, index=False)
    print(f"Wrote {args.output_table}", file=sys.stderr)

    # ── plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(max(10, len(samples) * 0.35), 7))

    cmap = plt.get_cmap("tab20")
    colors = {tag: cmap(i % 20) for i, tag in enumerate(pct.index)}

    bottom = pd.Series(0.0, index=samples)
    for tag in pct.index:
        values = pct.loc[tag]
        ax.bar(samples, values, bottom=bottom, label=tag, color=colors[tag])
        bottom = bottom + values

    ax.set_ylabel("Percent of reads (%)")
    ax.set_xlabel("Sample")
    title = "Contig/reference tag read-count composition per sample"
    if not args.include_unknown:
        title += " (excluding unknown)"
    ax.set_title(title)
    ax.set_ylim(0, 100)
    ax.tick_params(axis="x", rotation=90)
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8, title="tag")

    fig.tight_layout()
    fig.savefig(args.output)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
