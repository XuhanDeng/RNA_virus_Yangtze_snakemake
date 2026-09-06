# 03_RDRP_identification

## Purpose

Identifies RNA-dependent RNA polymerase (RdRP) candidates from assembled RNA virus scaffolds using a multi-tool strategy: HMM-based detection with RdRpCATCH, deep-learning classification with LucaProt (on ORFs), concatenation of all candidates, and motif A-D search using four complementary methods (HMM, PSI-BLAST, MMseqs2 profile search, and Diamond). Optionally, the same pipeline can be run on ICTV Riboviria reference sequences as a positive control.

## Rules

### Sample-level RdRP Detection

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `rdrpcatch` | HMM-based scan of scaffolds using all RdRpCATCH databases | Per-sample scaffolds FASTA | Annotated TSV, trimmed AA FASTA, full AA FASTA |
| `orffinder` | Predicts all ORFs from scaffolds using NCBI ORFfinder | Per-sample scaffolds FASTA | Per-sample ORF protein FASTA (`.faa`) |
| `lucaprot_rdrp` | Runs LucaProt deep-learning binary classifier on all ORFs | ORF FAA + `marker_db` | Per-sample CSV with RdRP predictions |
| `filter_lucaprot_threshold` *(optional)* | Re-filters LucaProt output to a stricter score threshold | LucaProt CSV (threshold 0.5) | LucaProt CSV at higher threshold |
| `extract_lucaprot_proteins_persample` | Extracts sequences of LucaProt-positive proteins into a FASTA | Filtered LucaProt CSV | Per-sample LucaProt protein FASTA |
| `extract_rdrpcatch_full_proteins` | Extracts full-length (not trimmed) proteins for RdRpCATCH hits | RdRpCATCH TSV + full AA FASTA | Per-sample RdRpCATCH full-protein FASTA |
| `cat_all_candidate_proteins` | Concatenates RdRpCATCH and LucaProt protein FASTAs across all samples | All per-sample FASTAs | `all_candidates.faa` |

### Motif Profile Building (RVMT Sequence_Library, motifs 1–4)

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `motif_add_consensus` | Adds consensus sequences to MSA files using `hhconsensus` | RVMT MSA directory | Marker `.consensus_done` |
| `motif_build_hmm` | Builds per-cluster HMM profiles with `hmmbuild`; concatenates into a single `.hmm` DB | MSA `.afa` files | `profiles.hmm` |
| `motif_msa_to_stockholm` | Converts `.afa` MSA files to Stockholm format for MMseqs2 | MSA directory | `.sto` files + marker `.sto_done` |
| `motif_build_mmseqs_profile` | Builds MMseqs2 profile databases from Stockholm files | Stockholm files | MMseqs2 profile DB + marker `.build_done` |

### Motif Searching

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `motif_search_hmmsearch` | Searches candidate proteins against HMM profiles | `profiles.hmm` + `all_candidates.faa` | `hmmsearch/mot.{motif}.tsv` |
| `motif_search_mmseqs` | Profile-vs-sequence search using MMseqs2, parallelised over 4 jobs | MMseqs2 profiles + `all_candidates.faa` | `mmseqs/mot.{motif}.tsv` |
| `motif_search_psiblast` | PSI-BLAST search using each MSA cluster as query | MSA files + `all_candidates.faa` | `psiblast/mot.{motif}.tsv` |
| `motif_search_diamond` | Diamond blastp with MSA sequences (gaps stripped) as query | MSA files + `all_candidates.faa` | `diamond/mot.{motif}.tsv` |
| `motif_merge_filter` | Merges all motif search results; classifies proteins into tiers (1a ABCD, 1b CABD-depermuted, 2 three-motif, 3 no-motif) | All motif TSVs + `all_candidates.faa` | Tier FASTAs, summary TSVs, suspicious-order TSV |

### ICTV Reference Pipeline (conditional on `config["ICTV"]["use"]`)

Mirrors all sample-level and motif rules for ICTV Riboviria sequences. Rules are named with the `_ICTV` suffix (e.g., `rdrpcatch_ICTV`, `lucaprot_ICTV`, `motif_merge_filter_ICTV`). Reuses HMM and MMseqs2 profiles built for the sample pipeline.

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `rna_samples` | Sample list |
| `seqkit["min_length"]` | Minimum scaffold length (locates input scaffold FASTAs) |
| `rna_reformated_scaffolds_dir` | Scaffolds directory from workflow 02 |
| `rdrp_catch["db_dir"]` / `["seq_type"]` / `["db_options"]` / `["threads"]` / `["runtime"]` | RdRpCATCH parameters |
| `orffinder["min_length"]` / `["strand"]` / `["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | ORFfinder parameters |
| `lucaprot["install_dir"]` / `["db_dir"]` / `["marker_db"]` | LucaProt installation and model paths |
| `lucaprot_rdrp["threshold"]` / `["filter_threshold"]` / `["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` / `["gpu_id"]` / `["gpu_devices"]` / `["gpus"]` / `["truncation_seq_length"]` / `["time_str"]` / `["step"]` / `["print_per_number"]` | LucaProt prediction parameters |
| `motif_search["prefix"]` / `["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` / `["parallel_module"]` / `["split_memory_limit"]` | Motif search parameters |
| `ICTV["use"]` | Boolean; enables ICTV reference pipeline when `true` |
| `ICTV["fasta"]` | ICTV Riboviria sequences FASTA |

## Dependencies / Tools Used

- `rdrpcatch` — HMM-based RdRP scan (conda: `rdrp_catch.yaml`)
- `ORFfinder` (NCBI) — ORF prediction (conda: `orffinder.yaml`)
- LucaProt `predict_many_samples.py` — deep-learning RdRP classification (conda: `lucaprot.yaml`)
- `seqkit grep` — protein extraction by ID (conda: `seqkit-spade.yaml`)
- `hhconsensus` (HH-suite) — consensus addition to MSAs (conda: `hhsuite.yaml`)
- `hmmbuild` / `hmmsearch` (HMMER) — HMM profile building and searching (conda: `hmmer.yaml`)
- `esl-reformat` (Easel/HMMER package) — MSA format conversion (conda: `mafft_hmmer_seqkit.yaml`)
- `mmseqs convertmsa` / `msa2profile` / `createdb` / `search` / `convertalis` (MMseqs2) — profile DB and search (conda: `mmseqs2.yaml`)
- `makeblastdb` / `psiblast` (BLAST+) — PSI-BLAST database and search (conda: `blast.yaml`)
- `diamond makedb` / `diamond blastp` — fast protein alignment (conda: `diamond.yaml`)
- GNU `parallel` — parallelises per-MSA steps; loaded as an HPC module
- Python scripts: `filter_lucaprot_threshold.py`, `extract_lucaprot_proteins_persample.py`, `motif_merge_filter.py`

## Notes

- Motifs 1–4 correspond to RdRP palm domain motifs A, B, C, D.
- LucaProt supports GPU acceleration; when `lucaprot_rdrp["gpu_id"] >= 0`, `CUDA_VISIBLE_DEVICES` is set and a GPU SLURM allocation is added via `slurm_extra`.
- The optional `filter_lucaprot_threshold` rule is only compiled into the DAG when `lucaprot_rdrp["filter_threshold"]` is set in config (Python `if` block at parse time).
- The ICTV branch is entirely conditional on `config["ICTV"]["use"]` being truthy.
- GNU Parallel is attempted via `module load`; if unavailable, rules fall back to sequential execution.
- MMseqs2 motif search splits threads across 4 parallel profile-vs-sequence jobs (`par_jobs = 4`).
- PSI-BLAST silently skips clusters with fewer than 2 sequences.
- `motif_merge_filter` classifies proteins into four tiers: Tier 1a (canonical A-B-C-D), Tier 1b (C-A-B-D permutation, depermuted), Tier 2 (three motifs), Tier 3 (no confirmed motif).
- RVMT Sequence_Library path is hardcoded: `database/RVMT/RVMT_Zenodo_V4/RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library`.
