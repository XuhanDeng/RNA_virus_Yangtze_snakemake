## Task: Post-process Spearman correlation results in workflow/05_esvirtu_correlation_analysis.smk

### Step 1 — Comment out FastSpar
- Comment out all FastSpar rules (`fastspar_prepare`, `fastspar_run`, `fastspar_bootstrap`, `fastspar_pvalues`, `fastspar_parse`) and their corresponding `rule all` target line.

### Step 2 — Add filtering rule after `rule spearman_all_vs_all`
Input files (already produced by `spearman_all_vs_all`):
- `result/05_esvirtu_correlation_analysis/2_all_vs_all/subspecies.all_vs_all.spearman.filtered.tsv`
- `result/05_esvirtu_correlation_analysis/2_all_vs_all/species.all_vs_all.spearman.filtered.tsv`

Apply the following to each file and save results under:
`result/05_esvirtu_correlation_analysis/3_spearman_analysis/{subspecies,species}/`

#### Pre-processing (apply before splitting)
1. Remove duplicate pairs — keep only rows where `seed <= target` (alphabetical), since A↔B pairs are symmetric.
2. Strip the `unknown::acc:` prefix from `seed` and `target` columns (e.g. `unknown::acc:AnQing_D_33_0000002790` → `AnQing_D_33_0000002790`).

#### Pattern 1 — known × known pairs
- Both `seed` and `target` are known viruses (label does not start with `unknown::` or `acc:`).
- Split by `r_group` threshold: one file per threshold (0.6, 0.7, 0.8, 0.9).
- Output: `known_known_pair/r{threshold}.tsv`

#### Pattern 2 — known × unknown pairs
- One side is a known virus, the other is unknown (label starts with `unknown::` or `acc:`).
- Split by `r_group` threshold: one file per threshold (0.6, 0.7, 0.8, 0.9).
- Add a `contig_name` column containing the cleaned contig ID of the unknown side.
- Output: `known_unknown_pair/r{threshold}.tsv`

#### Additional columns (both patterns)
Both the `spearman_all_vs_all` script and the filter output must include:
- `seed_nonzero_samples`: number of samples with non-zero RPKMF for the seed virus.
- `target_nonzero_samples`: number of samples with non-zero RPKMF for the target virus.

#### Resource config
Use `config["small_job"]` for the filter rule (lightweight Python job).
