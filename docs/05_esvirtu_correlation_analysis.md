# 05_esvirtu_correlation_analysis

## Purpose

Performs Spearman correlation analysis on viral abundance tables (RPKMF) produced by workflow 04, to find co-occurring viruses across samples. Results are split into known-known and known-unknown virus pairs at multiple correlation thresholds, then annotated with RdRP classification categories from workflow 03.

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `spearman_all_vs_all` | Computes pairwise Spearman correlations for all virus pairs at assembly, subspecies, and species levels; filters by minimum sample presence, p-value, and r thresholds | RPKMF tables (one per level) | `{level}.all_vs_all.spearman.filtered.tsv` |
| `filter_spearman_pairs` | Splits filtered correlations into known-known and known-unknown pairs at r thresholds 0.6 / 0.7 / 0.8 / 0.9 (subspecies and species levels only) | Per-level filtered Spearman TSV | `{level}/{pair_type}/r{threshold}.tsv` for each combination |
| `annotate_rdrp_category` | Joins known-unknown pair tables with `tier_summary.tsv` from workflow 03 and adds a `rdrp_category` column (tier1a / tier1b / tier2 / tier3 / not_in_list) | Known-unknown TSV + tier_summary.tsv | `{level}/known_unknown_pair/r{threshold}.annotated.tsv` |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `correlation["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Resources for Spearman all-vs-all |
| `correlation["min_samples"]` | Minimum number of samples in which a virus must be detected |
| `correlation["thresholds"]` | List of r thresholds to apply |
| `correlation["p_threshold"]` | Maximum p-value for a correlation to be retained |
| `correlation["chunk_size_all_vs_all"]` | Chunk size for the chunked pairwise computation |
| `scripts["spearman_all_vs_all"]` | Path to the all-vs-all Spearman script |
| `small_job["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Resources for filtering and annotation steps |

## Dependencies / Tools Used

- Python (pandas, scipy) — Spearman correlation calculation (conda: `python.yaml`)
- Python scripts: `spearman_all_vs_all.py`, `filter_spearman_pairs.py`, `annotate_rdrp_category.py`

## Notes

- Input RPKMF tables are taken from `result/04_modified_esvirtue_ribodetector/es/Merge` (the ribodetector variant of workflow 04).
- Three abundance levels are analysed: `assembly` (contig-level), `subspecies`, and `species`. Pair filtering is only applied at `subspecies` and `species` level.
- The four r thresholds (0.6, 0.7, 0.8, 0.9) are hardcoded as module-level constants `_THRESHOLDS`; only these exact values match the wildcard constraint `"0\\.6|0\\.7|0\\.8|0\\.9"`.
- A FastSpar (SparCC) pipeline is present in the file but entirely commented out. It was designed to use raw read-count tables rather than RPKMF; it can be re-enabled by uncommenting and adding its targets to `rule all`.
- `_TIER_SUMMARY_TSV` path is hardcoded: `result/03_RDRP_identification/4_motif_search/motif_results/tier_summary.tsv`. This is produced by the `motif_merge_filter` rule in workflow 03. The `rdrp_category` values correspond to motif-based tiers: **tier1a** (ABCD canonical), **tier1b** (CABD depermuted), **tier2** (3 motifs), **tier3** (tool-identified, <3 motifs or no motifs).
- Input RPKMF files use `ancient()` so stale timestamps do not trigger reruns.

## Run Command

```bash
nohup snakemake --snakefile workflow/05_esvirtu_correlation_analysis.smk \
    --executor slurm \
    --jobs 8 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/05_esvirtu_correlation_analysis/snakemake.log 2>&1 &
```
