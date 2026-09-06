# 101_kraken_check

## Purpose

Performs taxonomic classification of fastp-cleaned RNA reads using Kraken2, then re-estimates species-level abundances with Bracken at six taxonomic ranks (phylum through species). Per-sample Bracken reports are copied into a merged directory and then consolidated into cross-sample abundance tables by a Python script.

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `kraken2` | Classifies paired reads against the Kraken2 database at the configured confidence threshold; raw classification output is deleted after the report is written | fastp cleaned R1/R2 + Kraken2 DB | Per-sample Kraken2 report (`.report`) |
| `bracken_build` | Builds Bracken database files (k-mer distributions) from the Kraken2 DB at the configured k-mer and read lengths; creates a marker file on completion | Kraken2 DB directory | `bracken_build.done` marker |
| `bracken` | Runs Bracken six times per sample (ranks P, C, O, F, G, S) to re-estimate abundance at each rank | Kraken2 report + `bracken_build.done` | Six Bracken reports per sample (`_phylum`, `_classes`, `_orders`, `_family`, `_genus`, `_species`) |
| `merge_bracken` | Copies all per-sample Bracken reports into a shared directory organised by rank subdirectory (P/C/O/F/G/S) | All per-sample Bracken reports | Directory with rank subdirectories containing all sample reports |
| `merge_bracken_report` | Runs a Python script to consolidate per-rank report directories into cross-sample abundance tables | Merged report directory | Directory of merged tables (one per rank: `phylum`, `classes`, `orders`, `family`, `genus`, `species`) |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `rna_samples` | Sample list |
| `rna_fastp_dir` | Directory of fastp-cleaned reads (input to Kraken2) |
| `kraken2_dir` | Output directory for Kraken2 and Bracken per-sample reports |
| `kraken2_db["dir"]` | Path to the Kraken2 database (also used for Bracken) |
| `kraken2["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` / `["confidence"]` | Kraken2 parameters |
| `bracken_build["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` / `["kmer_length"]` / `["read_length"]` | Bracken database build parameters |
| `bracken["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` / `["read_length"]` | Bracken abundance estimation parameters |
| `merge_bracken["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Report copying parameters |
| `bracken_merge_dir` | Intermediate directory holding rank-organised Bracken reports |
| `bracken_report_dir` | Final output directory for consolidated cross-sample tables |
| `merge_bracken_report["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Merge script parameters |

## Dependencies / Tools Used

- `kraken2` — taxonomic classification (conda: `kraken2.yaml`)
- `bracken-build` — Bracken k-mer distribution construction (conda: `kraken2.yaml`)
- `bracken` — Bayesian abundance re-estimation (conda: `kraken2.yaml`)
- Python script `scripts/101_kraken_check/merge_profiling_reports.py` — cross-sample table consolidation (conda: `python.yaml`)

## Notes

- The `kraken2` rule writes a raw classification output file (`{output.out}`) and immediately deletes it with `rm -f`, keeping only the summary report to save disk space. Note: the rule references `{output.out}` but this output is not declared in the `output:` block — this is a latent bug; the `rm` will fail silently or target the wrong path.
- `bracken_build` only needs to run once per database; the `bracken_build.done` marker prevents redundant rebuilds.
- `merge_bracken` uses shell `for` loops to copy reports rather than Snakemake expand, so it has no per-file dependency tracking — any change to any report reruns the whole copy step.
- The `rule all` target is `directory(config["bracken_report_dir"])`, so Snakemake considers the workflow complete when that directory exists, regardless of its contents.
- `merge_bracken_report` uses `os.path.abspath(output.dir)` to construct the output path, ensuring the Python script receives an absolute path even if Snakemake changes working directory.
