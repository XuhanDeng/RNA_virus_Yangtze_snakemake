# 06_contig_bait

## Purpose

Recovers additional RNA virus contigs from the "no-RdRP" contig set by using known viral sequences as bait. RVMT NovoContigs and confirmed RdRP-positive contigs are used as queries in a nucleotide MMseqs2 search against contigs that passed no RdRP filter; hits passing stringent identity and coverage thresholds are extracted as a new recovered contig set.

## Rules

| Rule | Description | Key Inputs | Key Outputs |
|------|-------------|-----------|-------------|
| `combine_bait` | Concatenates RVMT NovoContigs and any-RdRP confirmed contigs into a single bait FASTA | `NovoContigs.fasta` + `any_rdrp.fasta` | `bait_combined.fasta` |
| `mmseqs_createdb_query` | Builds an MMseqs2 sequence database from the bait FASTA | `bait_combined.fasta` | `bait_db.index` |
| `mmseqs_createdb_target` | Builds an MMseqs2 sequence database from the no-RdRP contigs | `no_rdrp.fasta` | `target_db.index` |
| `mmseqs_search` | Nucleotide-vs-nucleotide MMseqs2 search (non-sensitive, -s 1) with loose initial filters: min alignment 120 bp, min identity 66%, target coverage 85% | Bait DB + Target DB | `result_db.index` |
| `mmseqs_convertalis` | Converts MMseqs2 result to a tabular TSV with alignment metrics | Bait DB + Target DB + Result DB | `hits.tsv` |
| `filter_bait_hits` | Applies stringent post-filters (e-value < configured threshold, identity > configured minimum, target coverage ≥ configured minimum); annotates each hit with its bait RdRP category; produces a best-bait-per-target map | `hits.tsv` + RdRP category lists | `hits.filtered.tsv` + `bait_target_map.tsv` |
| `extract_bait_contigs` | Extracts recovered contig sequences from `no_rdrp.fasta` using the filtered hit target IDs | `hits.filtered.tsv` + `no_rdrp.fasta` | `bait_recovered.txt` + `bait_recovered.fasta` |

## Key Config Parameters

| Config key | Description |
|-----------|-------------|
| `contig_bait["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | MMseqs2 search resources |
| `contig_bait["filter_evalue"]` | Maximum e-value for post-filter (e.g., 1e-9) |
| `contig_bait["filter_min_pident"]` | Minimum percent identity for post-filter (e.g., 95) |
| `contig_bait["filter_min_tcov"]` | Minimum target coverage for post-filter (e.g., 0.95) |
| `seqkit["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Resources for seqkit extraction step |
| `small_job["threads"]` / `["memory"]` / `["runtime"]` / `["partition"]` / `["account"]` | Resources for combine and filter steps |

## Dependencies / Tools Used

- `mmseqs createdb` / `search` / `convertalis` — nucleotide sequence database and search (conda: `mmseqs2.yaml`)
- `seqkit grep` — extract sequences by ID list (conda: `seqkit-spade.yaml`)
- Python script `filter_bait_hits.py` — post-filter and annotate hits (conda: `python.yaml`)

## Notes

- Bait sequences are sourced from two hardcoded paths:
  - `database/RVMT/RVMT_Zenodo_V4/RVMT_Zenodo_V4/Expanded_Contig_Set/NovoContigs.fasta`
  - `result/03_RDRP_identification/5_RdRp_contig/any_rdrp.fasta`
- Target sequences are sourced from: `result/03_RDRP_identification/7_no_RDRP_contig/no_rdrp.fasta`
- MMseqs2 is run in nucleotide mode (`--search-type 3`) at non-sensitive speed (`-s 1`) with `--cov-mode 1` (target coverage); this is a loose first pass — stringent filtering is done in `filter_bait_hits`.
- `filter_bait_hits.py` receives six RdRP category list files (High-confidence, RDRPCatch_palmscan, RDRPCatch_lucaprot, RDRPCatch_only, lucaprot_palmscan, lucaprot_only) to annotate each bait with its origin category.
- All bait and target inputs use `ancient()` to avoid reruns due to upstream timestamp changes.
- The `bait_target_map.tsv` output gives one row per recovered contig showing the best-matching bait and its RdRP category.
