# 07_RDRP_tree_construction

## Purpose

Builds a phylogenetic tree of RdRP proteins from Yangtze River contigs and ICTV reference sequences. The pipeline collects per-category protein FASTAs, de-permutes C-A-B-D permuted sequences to canonical A-B-C-D order, merges with RVMT reference sequences, aligns with MAFFT, trims with TrimAl, and infers a tree with IQ-TREE2. Steps 4–6 (alignment, trimming, tree) are currently commented out in `rule all` and the rules themselves are commented out.

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `collect_contig_input` | Copies per-category contig protein FASTAs into the workflow input directory | `result/03_RDRP_identification/8_RdRp_protein/{cat}.faa` | `result/07_phylogenetic_tree/0_input/Contig/{cat}.faa` |
| `collect_ictv_input` | Copies per-category ICTV protein FASTAs into the workflow input directory | `/scratch/.../8_RdRp_protein/{cat}.faa` | `result/07_phylogenetic_tree/0_input/ICTV/{cat}.faa` |
| `depermute` | De-permutes all Contig and ICTV FASTAs (C-A-B-D → A-B-C-D); merges all into a single combined FASTA and writes a summary TSV | Per-category Contig + ICTV FASTAs | Per-category depermuted FASTAs, `Contig_ICTV_rdrp.faa`, `Contig_ICTV_rdrp_summary.tsv` |
| `add_reference` | Concatenates depermuted Contig+ICTV proteins with RVMT RdRP reference sequences | `Contig_ICTV_rdrp.faa` + `RVMT_RdRp.faa` | `Contig_ICTV_RVMT_rdrp.faa` |
| `mafft_align` *(commented out)* | Multiple sequence alignment | Combined FASTA | Aligned FASTA |
| `trimal` *(commented out)* | Trims poorly aligned columns with automated mode | Aligned FASTA | Trimmed FASTA |
| `iqtree` *(commented out)* | Infers ML phylogenetic tree with model testing and 1000 ultrafast bootstraps | Trimmed FASTA | `.treefile` + support files |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `small_job["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Resources for collect, depermute, and add-reference rules |
| `mafft["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | MAFFT alignment resources (used when uncommented) |
| `iqtree["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | IQ-TREE2 resources (used when uncommented) |

## Dependencies / Tools Used

- Python script `02_depermute.py` — de-permutation of C-A-B-D sequences
- Python script `03_add_reference.py` — merging with RVMT references
- `mafft --auto` — multiple sequence alignment (conda: implied by `mafft` config block)
- `trimal -automated1` — alignment trimming
- `iqtree2 -m TEST -B 1000` — maximum-likelihood phylogenetic tree with model selection and ultrafast bootstrap

## Notes

- Four protein categories are processed: `High-confident`, `RDRPCatch_palm`, `lucaprot_palmscan`, `lucaprot_palmscan_extra`. These match subdirectories in `result/03_RDRP_identification/8_RdRp_protein/`.
- ICTV protein FASTAs are sourced from a **hardcoded absolute scratch path**: `/scratch/xddeng/yangtze/RNA/my_rna/result/03_RDRP_identification/ICTV/8_RdRp_protein`. This path will not work on a different user's account or cluster.
- RVMT RdRP reference path is also hardcoded: `database/RVMT/RVMT_Zenodo_V4/RVMT_Zenodo_V4/RVMT_RdRp.faa`.
- Steps 3–6 (`add_reference`, `mafft_align`, `trimal`, `iqtree`) are excluded from `rule all` (commented out). Only steps 1–2 (collect + depermute) are currently active targets.
- The `depermute` rule processes all categories in a single job rather than per-category, so any change to any category reruns all depermutation.
- Resource keys for `collect_contig_input` / `collect_ictv_input` use `mem_mb` (not `mem_mb_per_cpu`) — ensure the SLURM profile handles this key correctly.
