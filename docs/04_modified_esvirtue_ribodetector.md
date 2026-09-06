# 04_modified_esvirtue_ribodetector

## Purpose

Detects and quantifies RNA viruses per sample using a modified two-pass EsViritu strategy, where input reads are rRNA-depleted reads from **ribodetector** (rather than raw fastp output). Results are written under `result/04_modified_esvirtue_ribodetector/`. The workflow is structurally identical to `04_modified_esvirtue.smk` except for the read source and output paths.

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
| `esviritu_mapped_to_standard_database` | First pass: tags read IDs, runs EsViritu against ESViritu v3.2.4 standard DB; records mapped read IDs | ribodetector non-rRNA reads + ESViritu DB | Assembly summary, consensus, coverage, tagged reads, mapped read list |
| `extract_unmapped_reads_esviritu` | Extracts reads not mapped in first pass; deletes tagged reads on completion | Tagged reads + mapped read list | Unmapped R1/R2 FASTQ |
| `esviritu_mapped_to_assembly_database` | Second pass: runs EsViritu on unmapped reads against the custom contig database | Unmapped reads + custom DB (`.mmi` + metadata) | Assembly summary, consensus, coverage |
| `merge_esviritu_assembly_summary` | Concatenates first- and second-pass summaries per sample | First + second TSVs | Merged TSV |
| `merge_esviritu_rpkmf_all_samples` | Merges all sample TSVs; computes RPKMF, subspecies/species aggregations, and read counts | Per-sample merged TSVs | Multi-sample RPKMF TSV, subspecies/species TSVs, read-count tables |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `rna_samples` | Sample list |
| `rna_reformated_scaffolds_dir` | Scaffolds from workflow 02 |
| `ribodetector_dir` | Directory of ribodetector non-rRNA reads (input to first EsViritu pass) |
| `seqkit["min_length"]` | Minimum contig length for database |
| `merge_esviritu_database_1["min_length"]` / `["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Database merging |
| `databases["esviritu_db_v3.2.4"]` | Standard ESViritu database directory |
| `esviritu_map["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` / `["keep"]` / `["quiet"]` / `["extra"]` | EsViritu mapping parameters |
| `scripts["ref_cluster_select_assembly_only"]` / `["merge_es_tsv_assembly_only"]` / `["esviritu_merge_rpkmf"]` | Custom Python script paths |
| `ref_cluster_select["db_id_col"]` / `["sample_order"]` / `["sample_order_file"]` / `["with_cluster_id"]` | Optional cluster representative selection arguments |

## Dependencies / Tools Used

- `seqkit seq` / `seqkit grep` / `seqkit fx2tab` — sequence operations (conda: `seqkit-spade.yaml`)
- `makeblastdb` / `blastn` — clustering BLAST (conda: `checkv.yaml`)
- `minimap2 -d` — genome index (conda: `minimap2.yaml`)
- `EsViritu` — virus detection via read mapping (conda: `esviritu_map.yaml`)
- `samtools view` — extract mapped read IDs (conda: `esviritu_map.yaml`)
- `gawk` — read ID tagging
- Python scripts: `ref_cluster_select_assembly_only.py`, `merge_es_tsv_assembly_only.py`, `esviritu_merge_rpkmf.py`

## Notes

- Key difference from `04_modified_esvirtue.smk`: reads are taken from `config["ribodetector_dir"]` (rRNA-depleted) rather than `config["rna_fastp_dir"]` (fastp only).
- Key difference from `04_modified_esvirtue_use_nonrna.smk`: the first-pass EsViritu call uses the standard `EsViritu` CLI directly (not the custom Python wrapper `modified_esviritu`).
- This variant adds a `rpkmf.species.tsv` output to `rule all` that is absent from `04_modified_esvirtue.smk`.
- The `merge_esviritu_rpkmf_all_samples` rule here accepts raw merged TSVs (not pre-computed RPKMF) as input; RPKMF computation is done inside the merge script rather than as a separate step.
- `extract_unmapped_reads_esviritu` deletes tagged input reads after extraction to save disk.
- All output paths use the prefix `result/04_modified_esvirtue_ribodetector/` to distinguish from other 04 variants.
