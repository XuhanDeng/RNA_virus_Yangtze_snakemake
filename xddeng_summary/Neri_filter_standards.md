# Neri's Filter Standards for RdRP Motif Identification

Filters are applied in three successive stages. Each stage is stricter than the one before.

---

## Stage 1: Search-time cutoffs (inside runner scripts)

Applied during the search itself, before any results are read into R.
These are intentionally loose — the goal is to capture everything and let R do the real filtering.

| Tool | E-value | Query cover | Identity | Other |
|------|---------|------------|---------|-------|
| HHsearch | ≤ 0.5 | — | — | `-Z 100000` (return up to 100,000 hits) |
| HMMER hmmsearch | ≤ 0.5 | — | — | `--noali` |
| PSI-BLAST | ≤ 0.5 | ≥ 1% | — | `-word_size 2 -threshold 9 -dbsize 20000000` |
| MMseqs2 | ≤ 0.5 | — | — | `-k 6 -s 7.5 --max-seqs 6000` |
| Diamond BLASTp | ≤ 0.05 | ≥ 50% | ≥ 60% | `--more-sensitive -b3.0 -k 0` |

> Diamond is the only tool with a strict cutoff at search time (identity ≥ 60%, qcov ≥ 50%), because it searches singletons (short, specific sequences) rather than full profiles.

---

## Stage 2: R parsing filters (applied per search tool after loading)

Applied in `Annotations.r` immediately after `GenericHitsParser()` reads each hit table.
These are the **primary quality gates**.

### For RdRP hits (NrpsVSNVPCp block, `Annotations.r` lines 907–909)

```r
hits <- hits[evalue    <  1e-5 ]
hits <- hits[score     >  10   ]
hits <- hits[ali_Qcov_len > 10 ]   # alignment length on query > 10 AA
```

### For the full-sequence domain annotation databases (reference)

These thresholds vary slightly by database but follow the same pattern:

| Database | E-value | Score | Ali length | pCoverage |
|----------|---------|-------|-----------|-----------|
| RNAVirDB (HMM) | ≤ 1e-5 | ≥ 10 | ≥ 25 AA | — |
| Pfam34 (HMM) | ≤ 1e-6 | ≥ 10 | ≥ 25 AA | — |
| CDD (HMM) | < 1e-5 | ≥ 12 | — | — |
| ECOD (HMM) | < 1e-7 | ≥ 9 | — | — |
| SCOPe (HMM) | < 1e-6 | ≥ 9 | — | — |
| InterDomain (HMM) | < 1e-9 | ≥ 10 | — | ≥ 5% |
| CATH (HMM) | < 1e-6 | ≥ 8 | — | — |

> For motif identification specifically, the RNAVirDB thresholds (`evalue ≤ 1e-5`, `score ≥ 10`, `ali_len ≥ 25`) are the most directly relevant.

---

## Stage 3: Contig-level (presence/absence) filtration

Applied in `presence_absence_sys_args.r`. This stage removes contigs that match cellular organisms, used during the discovery pipeline to discard non-viral sequences. **This stage is about the contig, not the motif hit directly.**

### DNA filtration (remove cellular matches)

```
Diamond BLASTx vs NR:
  → keep contigs with NO hit to Bacteria/Archaea/Eukaryota
  → discard if: pident > 45%, alignment length > 30 AA

BLASTn vs NT:
  → discard cellular matches with bitscore > 82–85

DC-MegaBlast vs NT:
  → cull to best bitscore per contig, then discard cellular hits
```

Keyword exclusions applied before discarding (sequences mentioning these terms in their hit title are kept even if cellular taxonomy):
```
"virus", "phage", "capsid", "coat", "tail", "LysM", "lysis",
"viral", "vector", "AngRem", "PRSV", "ds-RNA element", "dsRNA element"
```

---

## Stage 4: "Tree RdRP" / RCR90 completeness criteria

These are the strictest filters, applied when selecting sequences for phylogenetic analysis. Documented in `Simplified_RdRPs/README.txt` and `Annotations.r`.

### Hit table level (`Set02_merged.tsv`)

```
E-value ≤ 1e-3     (one best hit per region per reading frame)
Culled to:  1 hit per region per reading frame
```

### Motif completeness (from `ite012_unfiltered.tsv`, implicit threshold)

Looking at the actual data range in `ite012_unfiltered.tsv` (267,731 rows):

```
score range:    11.6  →  highest values
evalue range:   3.95e-48  →  0.5 (unfiltered file keeps everything)
ali_len range:  ~20 AA minimum visible in data
pCoverage:      ~0.4 minimum visible in data
```

For sequences to qualify for the phylogenetic tree (RCR90 set), the README states criteria were **"far stricter"** than the disclosed hit table, requiring:

| Criterion | Requirement |
|-----------|------------|
| All four motifs present | Must have hits for mot.1 + mot.2 + mot.3 + mot.4 |
| Motif order identifiable | Positions q1 must be reliably assignable |
| Sequence completeness | Must cover enough of the RdRP domain core |
| Permutation resolvable | Permuted sequences must be de-permutable (C loop identifiable and cuttable) |
| No frameshifts preventing assembly | Must be de-frameshiftable to a single reading frame |

> Sequences that were too fragmented, too divergent, or where motif positions could not be reliably pinpointed were excluded from the tree set even if they passed the E-value threshold.

---

## Summary: recommended thresholds for your own analysis

If you are replicating Neri's approach on your candidate RdRPs:

### Minimum thresholds to accept a motif hit

```r
evalue    < 1e-5       # hard cutoff
score     > 10         # bit score
ali_len   > 10         # alignment length on your sequence (AA)
```

### Additional recommended filters

```r
pCoverage >= 0.40      # at least 40% of the motif profile is covered
                       # (implied by the data distribution in ite012_unfiltered.tsv)
```

### For inclusion in phylogenetic analysis (RCR90-equivalent)

```
1. All four motifs (mot.1–4) must have at least one passing hit
2. Motif order must be unambiguous (clear q1 ranking across the four motifs)
3. Total RdRP region length must be sufficient to cover the palm domain core
4. If permuted: the C-containing loop must be identifiable and excisable
5. Cull to 1 best hit per motif per sequence (highest score)
```

### What the final `ite012_unfiltered.tsv` represents

This file is **post Stage 2, pre Stage 4** — it has been quality-filtered by E-value/score/length but has NOT been filtered for completeness (all four motifs) or tree-inclusion criteria. It is the "disclosed for record" set, not the tree set. The tree set (RCR90 / `RdRPs_ali822x.faa`) is a strict subset where all four motifs were confirmed and the sequence was de-permuted/de-frameshifted.
