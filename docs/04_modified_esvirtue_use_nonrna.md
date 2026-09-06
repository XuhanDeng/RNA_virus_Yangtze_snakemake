# 04_modified_esvirtue_use_nonrna

## Purpose

Detects and quantifies RNA viruses per sample using a modified two-pass EsViritu strategy, where input reads are rRNA-depleted reads from **ribodetector** and the first EsViritu pass uses a **custom Python wrapper** (`modified_esviritu`) rather than the standard EsViritu CLI. Results are written under `result/04_modified_esvirtue/`.

## Rules

### Database Construction

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `merge_esviritu_database_1` | Concatenates all per-sample scaffolds (min-length filtered) into one merged FASTA | Per-sample renamed scaffolds | `esviritu_merged_db_only_assembly.fasta` |
| `merge_esviritu_database_2_checkv` | (Legacy/reference) All-vs-all BLASTn + ANI + CheckV-style clustering | Merged FASTA | BLAST DB, blast TSV, ANI TSV, cluster TSV |
| `select_cluster_representatives` | Selects one representative per cluster | Cluster TSV | `cluster_representatives.txt` |
| `extract_final_virus_pathogen_database` | Extracts representative sequences from merged FASTA | Representatives list + merged FASTA | `virus_pathogen_database.fna` |
| `index_final_virus_pathogen_database` | Builds a minimap2 index of the final database | `virus_pathogen_database.fna` | `virus_pathogen_database.mmi` |
| `seqkit_length_from_merged_db` | Computes per-sequence lengths | Merged FASTA | `length.txt` |
| `merge_esviritu_metadata` | Generates metadata TSV for representative sequences | Length file + IDs + FASTA | `virus_pathogen_database.all_metadata.tsv` |

### Virus Detection (Two-Pass EsViritu)

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `esviritu_mapped_to_standard_database` | First pass: tags read IDs, runs the **modified EsViritu Python script** against ESViritu v3.2.4 standard DB; records mapped read IDs | ribodetector non-rRNA reads + ESViritu DB | Assembly summary, consensus, coverage, tagged reads, mapped read list |
| `extract_unmapped_reads_esviritu` | Extracts reads not mapped in first pass; does **not** delete tagged reads (unlike other variants) | Tagged reads + mapped read list | Unmapped R1/R2 FASTQ |
| `esviritu_mapped_to_assembly_database` | Second pass: runs modified EsViritu on unmapped reads against the custom contig database (directory path, not `.mmi`) | Unmapped reads + custom DB directory | Assembly summary, consensus, coverage |
| `merge_esviritu_assembly_summary` | Concatenates first- and second-pass summaries per sample | First + second TSVs | Merged TSV |
| `rpkmf_esviritu_assembly_summary` | Normalises to RPKMF | Merged TSV | RPKMF TSV per sample |
| `merge_esviritu_rpkmf_all_samples` | Merges RPKMF across all samples; outputs subspecies/species and read-count tables; cleans up first_filter directory | Per-sample RPKMF TSVs | Multi-sample RPKMF, subspecies/species TSVs, read-count tables |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `rna_samples` | Sample list |
| `rna_reformated_scaffolds_dir` | Scaffolds from workflow 02 |
| `ribodetector_dir` | Directory of ribodetector non-rRNA reads |
| `seqkit["min_length"]` | Minimum contig length for database |
| `databases["esviritu_db_v3.2.4"]` | Standard ESViritu database directory |
| `scripts["modified_esviritu"]` | Path to the custom EsViritu Python wrapper |
| `scripts["ref_cluster_select_assembly_only"]` / `["merge_es_tsv_assembly_only"]` / `["esviritu_rpkmf"]` / `["esviritu_merge_rpkmf"]` | Custom Python script paths |
| `esviritu_map["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` / `["keep"]` / `["quiet"]` / `["extra"]` | EsViritu mapping parameters |
| `ref_cluster_select["db_id_col"]` / `["sample_order"]` / `["sample_order_file"]` / `["with_cluster_id"]` | Optional cluster representative selection arguments |
| `rpkmf_esviritu_assembly_summary["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | RPKMF normalisation parameters |
| `merge_esviritu_rpkmf_all_samples["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Cross-sample merge parameters |

## Dependencies / Tools Used

- `seqkit seq` / `seqkit grep` / `seqkit fx2tab` — sequence operations (conda: `seqkit.yaml`)
- `makeblastdb` / `blastn` — clustering BLAST (conda: `checkv.yaml`)
- `minimap2 -d` — genome index (conda: `minimap2.yaml`)
- Custom Python script (`modified_esviritu`) wrapping EsViritu logic (conda: `esviritu_map.yaml`)
- `samtools view` — extract mapped read IDs (conda: `esviritu_map.yaml`)
- `gawk` — read ID tagging
- Python scripts: `ref_cluster_select_assembly_only.py`, `merge_es_tsv_assembly_only.py`, `esviritu_rpkmf.py`, `esviritu_merge_rpkmf.py`

## Notes

- Distinguishing feature: both EsViritu passes use `python {params.script}` (the custom wrapper) rather than calling `EsViritu` directly. The custom script is configured via `config["scripts"]["modified_esviritu"]`.
- The second-pass database input is a **directory path** (`result/04_modified_esvirtue/databases/final_merged_database_only_assembly`) rather than an explicit `.mmi` file, unlike the ribodetector variant which passes `$(dirname {input.db_mmi})`.
- `extract_unmapped_reads_esviritu` does **not** delete input tagged reads in this variant (no `rm -f {input.r1} {input.r2}` in the shell command).
- `merge_esviritu_rpkmf_all_samples` deletes `result/04_modified_esvirtue/es/first_filter` on completion.
- Uses `../envs/seqkit.yaml` (not `seqkit-spade.yaml`) for seqkit steps.
- Output paths use `result/04_modified_esvirtue/` — the same prefix as `04_modified_esvirtue.smk`. Running both workflows would cause output conflicts; they are mutually exclusive alternatives.
