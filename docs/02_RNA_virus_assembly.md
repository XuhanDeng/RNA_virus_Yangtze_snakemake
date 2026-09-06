# 02_RNA_virus_assembly

## Purpose

Processes paired-end RNA-seq reads for each sample through a four-step pipeline: quality control (fastp), rRNA removal (ribodetector), metagenomic assembly (SPAdes), and contig renaming/length filtering. Produces per-sample scaffolds with standardised sequence IDs ready for downstream virus identification.

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `fastp_qc` | Adapter trimming, poly-G/X trimming, and quality filtering; retains both paired and unpaired reads | `{sample}.R1/R2.fq.gz` | Paired cleaned reads (`_1P/_2P`), unpaired reads (`_U1/_U2`), HTML/JSON QC reports |
| `ribodetector_rrna_removal` | Removes ribosomal RNA reads using a deep-learning CPU classifier | Paired cleaned reads | Non-rRNA paired reads (`_nonrrna.1/2.fq.gz`) |
| `spades_assembly` | Metagenomic de-novo assembly with custom k-mer list; error correction disabled | Non-rRNA reads | `scaffolds.fasta` |
| `rename_filter_assemblies` | Filters scaffolds by minimum length and renames them `{sample}_{nr}` | `scaffolds.fasta` | `{sample}_scaffolds_rename_{min_len}.fasta` |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `rna_samples` | List of sample names |
| `rna_input_dir` | Directory containing raw RNA FASTQ files |
| `rna_fastp_dir` | Output directory for fastp results |
| `rna_fastp["quality_threshold"]` / `["length_required"]` / `["threads"]` / `["memory_per_cpu"]` / `["runtime"]` / `["partition"]` / `["account"]` | fastp parameters |
| `ribodetector_dir` | Output directory for ribodetector results |
| `ribodetector["min_length"]` / `["chunk_size"]` / `["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | ribodetector_cpu parameters |
| `rna_spades_dir` | Output directory for SPAdes assemblies |
| `rna_spades["k_values"]` / `["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | SPAdes parameters |
| `rna_reformated_scaffolds_dir` | Output directory for renamed/filtered scaffolds |
| `seqkit["min_length"]` / `["nr_width"]` / `["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Contig length filter and renaming parameters |

## Dependencies / Tools Used

- `fastp` — quality control (conda: `fastp.yaml`)
- `ribodetector_cpu` — rRNA read removal (conda: `ribodetector.yaml`)
- `spades.py` — metagenomic assembly with `--meta --only-assembler` (conda: `seqkit-spade.yaml`)
- `seqkit seq` / `seqkit replace` — length filtering and sequence renaming (conda: `seqkit-spade.yaml`)

## Notes

- Raw input files use `ancient()` so stale timestamps on raw data do not trigger reruns.
- SPAdes memory is computed dynamically from `mem_mb_per_cpu * cpus_per_task // 1024` (MB → GB conversion).
- The assembly output path includes the suffix `_no_correction` because error correction is disabled (`--only-assembler`).
- Final scaffold IDs follow the pattern `{sample}_{nr}` where `{nr}` is zero-padded to `nr_width` digits.
- The `rule all` target is parameterised by `config["seqkit"]["min_length"]`, so all per-sample FASTAs use the same minimum length threshold simultaneously.
- Compare with `02_RNA_virus_assembly_sortmerna.smk` which substitutes ribodetector with SortMeRNA + BBTools repair.
