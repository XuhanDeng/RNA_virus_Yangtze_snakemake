# Bait contig with RDRP

method reference: Comprehensive identification of RNA virus contigs across metatranscriptomes

> Because metatranscriptome assemblies can often yield incomplete genomes that would not fulfil
> the criteria for de novo RdRP detection, we used the "VR1507" contig set to initiate a secondary
> "sweeping" scan for additional RNA virus contigs from the non-clustered, non-filtered bulk-set of
> metatranscriptomic contigs. To this end, "VR1507" was used as bait for highly similar contigs in
> the "bulk-set", using a non-sensitive mmseqs search
> (mmseqs search –search-type 3 –min-aln-len 120 –min-seq-id 0.66 -s 1 -c 0.85 –cov-mode 1)
> followed by stringent filtering of the recovered matches
> (E-value < 1e-9, Identity > 95%, target-Coverage ≥ 95%).
> The envelopment criterion (target coverage ≥ 95%) was added to avoid capture of chimeric or
> otherwise uncertain nucleic regions that extend beyond the query.

---

## Inputs

| Role | Path |
|------|------|
| Bait 1 — RVMT NovoContigs | `database/RVMT/RVMT_Zenodo_V4/RVMT_Zenodo_V4/Expanded_Contig_Set/NovoContigs.fasta` |
| Bait 2 — RdRP-confirmed contigs | `result/03_RDRP_identification/5_RdRp_contig/any_rdrp.fasta` |
| Target — no-RdRP contigs | `result/03_RDRP_identification/7_no_RDRP_contig/no_rdrp.fasta` |
| RdRP category lists | `result/03_RDRP_identification/4_merged/{High-confidence_RdRp,RDRPCatch_palmscan,RDRPCatch_lucaprot,RDRPCatch_only,lucaprot}.txt` |

## Output structure

```
result/06_contig_bait/
  1_bait_db/
    bait_combined.fasta          — RVMT + any_rdrp merged bait sequences
    bait_db*                     — mmseqs2 query database
  2_target_db/
    target_db*                   — mmseqs2 target database (no_rdrp.fasta)
  3_mmseqs_search/
    result_db*                   — raw mmseqs2 result database
    hits.tsv                     — all hits (query, target, pident, alnlen, evalue, tcov, ...)
  4_filtered_hits/
    hits.filtered.tsv            — hits passing: evalue<1e-9, pident>95%, tcov≥95%
    bait_target_map.tsv          — one row per (recovered_contig × bait_source):
                                   recovered_contig | bait_source | bait_category |
                                   best_bait_contig | pident | tcov | evalue | n_bait_hits
  5_bait_contigs/
    bait_recovered.txt           — unique recovered contig IDs
    bait_recovered.fasta         — their sequences extracted from no_rdrp.fasta
```

### bait_target_map.tsv columns

| column | description |
|--------|-------------|
| `recovered_contig` | target contig ID (from no_rdrp.fasta) |
| `bait_source` | `own_rdrp` (from any_rdrp.fasta) or `RVMT` (from NovoContigs) |
| `bait_category` | RdRP category of the best bait contig: `High-confidence_RdRp`, `RDRPCatch_palmscan`, `RDRPCatch_lucaprot`, `RDRPCatch_only`, `lucaprot`, or `RVMT` |
| `best_bait_contig` | query ID with the lowest e-value for this target × source pair |
| `pident` | percent identity of that best hit |
| `tcov` | target coverage of that best hit (0–1) |
| `evalue` | e-value of that best hit |
| `n_bait_hits` | total distinct bait contigs hitting this target (all sources combined) |

## Workflow file

`workflow/06_contig_bait.smk`

### Rules

| step | rule | tool |
|------|------|------|
| 1 | `combine_bait` | `cat` |
| 2a | `mmseqs_createdb_query` | mmseqs2 |
| 2b | `mmseqs_createdb_target` | mmseqs2 |
| 3 | `mmseqs_search` | mmseqs2 |
| 4 | `mmseqs_convertalis` | mmseqs2 |
| 5 | `filter_bait_hits` | `scripts/06_contig_bait/filter_bait_hits.py` |
| 6 | `extract_bait_contigs` | seqkit |

### mmseqs search parameters (from paper)

| parameter | value | meaning |
|-----------|-------|---------|
| `--search-type` | 3 | nucleotide–nucleotide |
| `--min-aln-len` | 120 | minimum alignment length |
| `--min-seq-id` | 0.66 | loose initial identity filter |
| `-s` | 1 | non-sensitive mode |
| `-c` | 0.85 | coverage threshold |
| `--cov-mode` | 1 | target coverage |

### Post-filter thresholds

| filter | value |
|--------|-------|
| E-value | < 1e-9 |
| Identity | > 95% |
| Target coverage | ≥ 95% |

## Scripts

- `scripts/06_contig_bait/filter_bait_hits.py`
  - Reads mmseqs2 `hits.tsv` (no header, 15 columns)
  - Applies stringent post-filters
  - Writes `hits.filtered.tsv` (all passing alignments)
  - Writes `bait_target_map.tsv` (best hit per target × bait_source, with RdRP category annotation)
  - Labels bait source by checking query ID against `any_rdrp.fasta` IDs
  - Labels bait category by checking query ID against the 5 RdRP category `.txt` files
