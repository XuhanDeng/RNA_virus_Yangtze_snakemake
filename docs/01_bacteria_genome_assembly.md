# 01_bacteria_genome_assembly

## Purpose

Assembles bacterial genomes from paired-end raw FASTQ reads for each sample: quality-trims with fastp, assembles with SPAdes (meta mode), renames contigs to a standardised `Bac_{sample}_{nr}` format, and finally clusters all samples together with MMseqs2 `easy-linclust` to remove redundant contigs.

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `fastp_quality_control` | Adapter trimming, poly-G/X trimming, and quality filtering of raw paired reads | `{sample}.R1/R2.raw.fastq.gz` | Cleaned reads (`_1P/_2P.fq.gz`), HTML/JSON QC reports |
| `spades_assembly_and_reformat` | Runs SPAdes in metagenomic mode, filters contigs by minimum length, renames them `Bac_{sample}_{10-digit-nr}`, and writes a renaming TSV | Cleaned reads | `{sample}_reformated.fa`, `{sample}_renaming.tsv`, `{sample}_original_contigs.fasta` |
| `cleanup_fastp` | Deletes cleaned FASTQ files after assembly is confirmed complete | `{sample}_reformated.fa` | `.cleaned` marker file |
| `seqkit_overXbp` | Alternative length-filter + renaming rule producing a per-sample FASTA of contigs above a configurable size threshold | `{sample}_original_contigs.fasta` | `{sample}_over{N}bp.fa` |
| `cat_all_scaffolds` | Concatenates all per-sample over-threshold FASTA files into one combined file | Per-sample over-Xbp FASTAs | `all_samples_over100bp.fa` |
| `mmseqs_easy_linclust` | Clusters all assembled scaffolds across samples to remove near-identical sequences | `all_samples_over100bp.fa` | `linclust_rep_seq.fasta` (cluster representatives) |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `bacteria_samples` | List of sample names to process |
| `bacteria_input_dir` | Directory containing raw input FASTQ files |
| `bacteria_fastp_dir` | Output directory for fastp results |
| `bacteria_fastp["threads"]` / `["qualified_quality_phred"]` / `["length_required"]` / `["runtime"]` | fastp run parameters |
| `bacteria_spades_dir` / `bacteria_reformated_scaffolds_dir` | SPAdes working and final output directories |
| `bacteria_spades["threads"]` / `["memory"]` / `["kmers"]` / `["memory_per_cpu"]` / `["runtime"]` / `["partition"]` | SPAdes resource and assembly parameters |
| `reformat_scaffolds["min_length"]` / `["nr_width"]` / `["runtime"]` | Minimum contig length and zero-padding width for renaming |
| `seqkit_overXbp["min_length"]` / `["nr_width"]` / `["threads"]` / `["runtime"]` / `["partition"]` | Parameters for the alternative length-filter rule |
| `easy-linclust_dir` | Output directory for MMseqs2 linclust results |
| `easy-linclust["min_seq_id"]` / `["cov_mode"]` / `["cluster_mode"]` / `["coverage"]` / `["threads"]` / `["memory_per_cpu"]` / `["runtime"]` / `["partition"]` | MMseqs2 clustering parameters |
| `regular_memory` / `regular_partition` / `account` | Default SLURM resource settings |
| `small_job["memory"]` / `["runtime"]` | Resource settings for the cleanup rule |

## Dependencies / Tools Used

- `fastp` — quality control and adapter trimming (conda: `fastp.yaml`)
- `spades.py` — metagenomic assembly with `--meta --only-assembler` (conda: `seqkit-spade.yaml`)
- `seqkit seq` / `seqkit replace` / `seqkit fx2tab` — sequence filtering and renaming (conda: `seqkit-spade.yaml`)
- `mmseqs easy-linclust` — sequence clustering to reduce redundancy (conda: `mmseqs.yaml`)

## Notes

- Raw input files are declared with `ancient()` so Snakemake does not rerun the pipeline if raw files are older than intermediate outputs.
- `spades_assembly_and_reformat` removes the SPAdes working directory (`{sample}` folder under `bacteria_spades_dir`) after the assembly is reformatted, saving disk space.
- `cleanup_fastp` removes cleaned reads after assembly to save storage; the `.cleaned` marker ensures this only happens once per sample.
- `seqkit_overXbp` and `cat_all_scaffolds` / `mmseqs_easy_linclust` are **not** in the default `rule all` targets — they must be requested explicitly or added to `rule all` for the clustering branch to run.
- Contig IDs are zero-padded to `nr_width` digits (e.g., `Bac_SampleA_0000000001`) to guarantee lexicographic sort order.
