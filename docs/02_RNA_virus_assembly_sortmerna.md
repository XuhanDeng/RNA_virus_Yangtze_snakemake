# 02_RNA_virus_assembly_sortmerna

## Purpose

Alternative RNA virus assembly workflow that uses SortMeRNA (database-based) instead of ribodetector (deep-learning) for rRNA removal. After SortMeRNA filtering, BBTools `repair.sh` is used to re-synchronise read pairs before SPAdes assembly. Final scaffold length is hard-coded to 1000 bp in the output paths.

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `fastp_qc` | Adapter trimming, poly-G/X trimming, quality filtering; retains paired and unpaired reads | `{sample}.R1/R2.fq.gz` | Paired cleaned reads (`_1P/_2P`), unpaired reads (`_U1/_U2`), HTML/JSON QC reports |
| `sortmerna_rrna_removal` | Maps reads against multiple rRNA reference databases and writes non-rRNA reads | Paired cleaned reads, multiple `--ref` DB FASTAs | Non-rRNA forward/reverse reads (`_nonrrna_fwd/_rev.fq.gz`) |
| `bbtools_repair` | Re-pairs and synchronises forward/reverse reads that may have become desynchronised after SortMeRNA | Non-rRNA reads | Re-paired reads (`_paired_1/_2.fq.gz`) |
| `spades_assembly` | Metagenomic de-novo assembly; error correction disabled | Re-paired reads | `scaffolds.fasta` |
| `rename_filter_assemblies` | Filters scaffolds to ≥ `seqkit.min_length` bp and renames them `{sample}_{nr}` | `scaffolds.fasta` | `{sample}_scaffolds_rename_1000.fasta` |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `rna_samples` | List of sample names |
| `rna_input_dir` | Directory containing raw RNA FASTQ files |
| `rna_fastp_dir` | Output directory for fastp results |
| `rna_fastp["quality_threshold"]` / `["length_required"]` / `["threads"]` | fastp parameters |
| `sortmerna_dir` | Output directory for SortMeRNA non-rRNA reads |
| `sortmerna["dbs"]` | List of rRNA reference FASTA paths (passed as multiple `--ref` flags) |
| `sortmerna["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | SortMeRNA resource parameters |
| `bbtools_dir` | Output directory for BBTools repair results |
| `bbtools_repair["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | BBTools resource parameters |
| `rna_spades_dir` | Output directory for SPAdes assemblies |
| `rna_spades["k_values"]` / `["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | SPAdes parameters |
| `rna_reformated_scaffolds_dir` | Output directory for renamed/filtered scaffolds |
| `seqkit["min_length"]` / `["nr_width"]` / `["threads"]` | Contig renaming parameters |

## Dependencies / Tools Used

- `fastp` — quality control (conda: `fastp.yaml`)
- `sortmerna` — rRNA read removal against multiple databases (conda: `sortmerna.yaml`)
- `repair.sh` (BBTools) — read pair resynchronisation (conda: `bbtools.yaml`)
- `spades.py` — metagenomic assembly with `--meta --only-assembler` (conda: `seqkit-spade.yaml`)
- `seqkit seq` / `seqkit replace` — length filtering and renaming (conda: `seqkit-spade.yaml`)

## Notes

- There is a hardcoded `--threads 10` in the `sortmerna_rrna_removal` shell command along with a comment `# remeber change --threads 10 to {threads} after testing` — the thread count is not actually taken from the config for SortMeRNA.
- SortMeRNA is run with `--fastx --out2 --paired-in` flags, outputting separate forward/reverse non-rRNA FASTQ files; a per-sample `workdir` is used for its internal index.
- The output scaffold FASTA path contains the literal string `rename_1000` (not a wildcard), because the minimum length is fixed at 1000 bp in this variant.
- The `rule all` target only covers the renamed scaffolds at 1000 bp; the BBTools repair and SortMeRNA output files are intermediate and not explicitly requested.
- The `account` for fastp uses the global `config["account"]` (not `config["rna_fastp"]["account"]`) — a minor inconsistency with the ribodetector variant.
- Compare with `02_RNA_virus_assembly.smk` which uses ribodetector instead of SortMeRNA+BBTools.
