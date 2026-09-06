# 100_figure

## Purpose

Generates network visualisation figures of virus co-occurrence (Spearman correlation) results from workflow 05. For each combination of taxonomic level, correlation threshold, and display variant, it produces four output files: SVG with labels, SVG without labels, interactive HTML with labels, and interactive HTML without labels.

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `plot_network` | Runs a Python network plotting script that reads known-unknown (KU) and known-known (KK) correlation pair tables and renders a network graph at the specified threshold and display variant | KU annotated pair directory + KK pair directory | `ku{threshold}_with_label.svg`, `ku{threshold}_no_label.svg`, `ku{threshold}_with_label.html`, `ku{threshold}_no_label.html` |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `network_figure["levels"]` | List of taxonomic levels to plot (e.g., `subspecies`, `species`) |
| `network_figure["thresholds"]` | List of r thresholds to plot (e.g., `0.6`, `0.7`, `0.8`, `0.9`) |
| `network_figure["variants"]` | Dict of display variants; each variant has keys `kk_threshold`, `hide_not_in_list`, and `kk_only` |
| `small_job["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | SLURM resources for each plot job |

## Dependencies / Tools Used

- Python script `scripts/100_network_figure/plot_network.py` — network graph generation (conda: `python.yaml`; likely uses networkx, matplotlib or pyvis for HTML output)

## Notes

- Input directories use `ancient()` to avoid reruns when upstream correlation files are older than figure outputs.
- The `variant` wildcard drives three conditional CLI flags passed to the plot script:
  - `kk_threshold`: if `"all"`, passes `--all-kk`; if a numeric string, passes `--kk-threshold {value}`; if empty, no KK threshold flag is added.
  - `hide_not_in_list`: if `"True"`, passes `--hide-not-in-list`.
  - `kk_only`: if `"True"`, passes `--kk-only`.
- Wildcard constraints enforce that `level` is `subspecies|species` and `threshold` matches `0.6|0.7|0.8|0.9`.
- KU inputs come from `result/05_esvirtu_correlation_analysis/4_spearman_analysis_RdRp_annotated/{level}/known_unknown_pair/`.
- KK inputs come from `result/05_esvirtu_correlation_analysis/3_spearman_analysis/{level}/known_known_pair/`.
- Each `plot_network` invocation reads all threshold files within the KU/KK directories; the `--min-threshold` argument selects which threshold to display.
