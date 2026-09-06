# 04_modified_esvirtue

## Purpose

Detects and quantifies RNA viruses per sample using a modified two-pass EsViritu strategy. Reads first pass through the standard ESViritu database; unmapped reads are then screened against a custom database built from assembled contigs. Results from both passes are merged, normalised to RPKMF, and aggregated across all samples. Input reads come from fastp output (not rRNA-depleted).

## Rules

### Database Construction

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `merge_esviritu_database_1` | Concatenates all per-sample scaffolds (filtered by minimum length) into one merged FASTA | Per-sample renamed scaffolds | `esviritu_merged_db_only_assembly.fasta` |
| `merge_esviritu_database_2_checkv` | (Legacy/reference) All-vs-all BLASTn + ANI calculation + CheckV-style clustering of merged DB | Merged FASTA | BLAST DB, blast TSV, ANI TSV, cluster TSV |
| `select_cluster_representatives` | Selects one representative sequence per cluster | Cluster TSV | `cluster_representatives.txt` |
| `extract_final_virus_pathogen_database` | Extracts representative sequences from merged FASTA | Representatives list + merged FASTA | `virus_pathogen_database.fna` |
| `index_final_virus_pathogen_database` | Builds a minimap2 index of the final database | `virus_pathogen_database.fna` | `virus_pathogen_database.mmi` |
| `seqkit_length_from_merged_db` | Computes sequence lengths from the merged database | Merged FASTA | `length.txt` |
| `merge_esviritu_metadata` | Generates a metadata TSV for all representative sequences | Length file + IDs + FASTA | `virus_pathogen_database.all_metadata.tsv` |

### Virus Detection (Two-Pass EsViritu)

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `esviritu_mapped_to_standard_database` | First pass: runs a modified EsViritu script against the standard ESViritu v3.2.4 database; tags read IDs for pair tracking; saves list of mapped read IDs | fastp cleaned reads + ESViritu DB dir | Assembly summary, consensus, coverage, tagged reads, mapped read list |
| `extract_unmapped_reads_esviritu` | Extracts reads that did NOT map in the first pass, using seqkit grep exclusion | Tagged reads + mapped read list | Unmapped R1/R2 FASTQ |
| `esviritu_mapped_to_assembly_database` | Second pass: runs modified EsViritu on unmapped reads against the custom contig database | Unmapped reads + custom DB (`.mmi` + metadata) | Assembly summary, consensus, coverage |
| `merge_esviritu_assembly_summary` | Concatenates first- and second-pass assembly summaries per sample | First + second TSVs | Merged TSV (with copies of each pass) |
| `rpkmf_esviritu_assembly_summary` | Normalises merged counts to RPKMF (reads per kilobase per million fragments) | Merged TSV | RPKMF TSV |
| `merge_esviritu_rpkmf_all_samples` | Merges RPKMF tables across all samples; also outputs subspecies/species aggregations and raw read count tables | Per-sample RPKMF TSVs | Multi-sample RPKMF, subspecies RPKMF, read-count tables |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `rna_samples` | Sample list |
| `rna_reformated_scaffolds_dir` | Scaffolds from workflow 02 |
| `seqkit["min_length"]` | Minimum contig length for database construction |
| `merge_esviritu_database_1["min_length"]` / `["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Database merging parameters |
| `checkv["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` / `["blast_outfmt"]` / `["blast_max_target_seqs"]` / `["min_ani"]` / `["min_coverage"]` / `["min_qcov"]` | CheckV clustering parameters (legacy rule) |
| `scripts["ref_cluster_select_assembly_only"]` / `["merge_es_tsv_assembly_only"]` / `["modified_esviritu"]` / `["esviritu_rpkmf"]` / `["esviritu_merge_rpkmf"]` | Custom Python script paths |
| `databases["esviritu_db_v3.2.4"]` | Path to standard ESViritu database directory |
| `esviritu_map["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` / `["keep"]` / `["quiet"]` / `["extra"]` | EsViritu mapping parameters |
| `ref_cluster_select["db_id_col"]` / `["sample_order"]` / `["sample_order_file"]` / `["with_cluster_id"]` | Optional cluster representative selection arguments |
| `index_final_virus_pathogen_database["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Minimap2 index parameters |
| `merge_esviritu_assembly_summary["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Merge parameters |
| `rpkmf_esviritu_assembly_summary["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | RPKMF normalisation parameters |
| `merge_esviritu_rpkmf_all_samples["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Cross-sample merge parameters |

## Dependencies / Tools Used

- `seqkit seq` / `seqkit grep` / `seqkit fx2tab` — sequence filtering and length calculation (conda: `seqkit-spade.yaml`)
- `makeblastdb` / `blastn` — all-vs-all nucleotide BLAST for clustering (conda: `checkv.yaml`)
- `minimap2 -d` — genome database indexing (conda: `minimap2.yaml`)
- `EsViritu` (via Python script `modified_esviritu`) — read mapping and virus detection (conda: `esviritu_map.yaml`)
- `samtools view` — extract mapped read IDs from BAM (conda: `esviritu_map.yaml`)
- `gawk` — read ID tagging for pair disambiguation
- Python scripts: `ref_cluster_select_assembly_only.py`, `merge_es_tsv_assembly_only.py`, `esviritu_rpkmf.py`, `esviritu_merge_rpkmf.py`

## Notes

- This variant takes reads from `rna_fastp_dir` (fastp output, **not** rRNA-depleted). Compare with `04_modified_esvirtue_ribodetector.smk` and `04_modified_esvirtue_use_nonrna.smk` which use ribodetector-filtered reads.
- The `merge_esviritu_database_2_checkv` (clustering) rule is kept for reference but is **not** in the default `rule all` targets; the workflow in practice uses only `merge_esviritu_database_1` + `select_cluster_representatives`.
- Read IDs are tagged with `_1` / `_2` suffixes before the first EsViritu pass to allow the unmapped-read extraction step to match read pairs correctly.
- `extract_unmapped_reads_esviritu` deletes the tagged FASTQ inputs after extraction to save disk space.
- `esviritu_mapped_to_assembly_database` deletes R1/R2 unmapped reads and the first-pass temp directory on completion.
- The `merge_esviritu_rpkmf_all_samples` rule deletes `result/04_modified_esvirtue/es/first_filter` on completion.
- The `build_ref_cluster_args` Python helper function at workflow top constructs optional CLI flags for the cluster representative selection script based on config keys.
