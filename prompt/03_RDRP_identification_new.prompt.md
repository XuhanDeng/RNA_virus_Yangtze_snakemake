# 03_RDRP_identification pipeline

## Overview

Four steps to identify RdRp candidates from RNA virus assemblies.
palm_annot has been removed — replaced by RVMT motif search (Step 4).

---

## Step 1 — RdRpCATCH (per sample + ICTV)

**Input:**
```
result/02_RNA_virus_assembly/4_rename_assembly/rename_1000/{sample}_scaffolds_rename_1000.fasta
```

**Output:**
```
result/03_RDRP_identification/1_RdRpCATCH/{sample}/{sample}_scaffolds_rename_1000_rdrpcatch_output_annotated.tsv
result/03_RDRP_identification/1_RdRpCATCH/{sample}/{sample}_scaffolds_rename_1000_rdrpcatch_fasta/
    {sample}_scaffolds_rename_1000_trimmed_aminoacid_sequences.fasta   ← trimmed RdRp AA (per sample)
    {sample}_scaffolds_rename_1000_full_aminoacid_sequences.fasta      ← full-length ORF AA (per sample)
```

ICTV output (scratch):
```
/scratch/xddeng/yangtze/RNA/my_rna/result/03_RDRP_identification/1_RdRpCATCH/ICTV/
    Riboviria_sequences_rdrpcatch_output_annotated.tsv
    Riboviria_sequences_rdrpcatch_fasta/Riboviria_sequences_trimmed_aminoacid_sequences.fasta
```

---

## Step 2 — ORFfinder + LucaProt (per sample + ICTV)

### Step 2a — ORFfinder

**Input:**
```
result/02_RNA_virus_assembly/4_rename_assembly/rename_1000/{sample}_scaffolds_rename_1000.fasta
```

**Output:**
```
result/03_RDRP_identification/3_lucaprot/{sample}/{sample}_orfs.faa
```

### Step 2b — LucaProt RdRp prediction

**Input:**
```
result/03_RDRP_identification/3_lucaprot/{sample}/{sample}_orfs.faa
```

**Output:**
```
result/03_RDRP_identification/3_lucaprot/{sample}/RdRPs_only_using_threshold0.900000.csv
```

Columns: `protein_id, seq, prob, label`
ID format: `>lcl|ORF{n}_{sample}_{contig}:{start}:{end} unnamed protein product`

### Step 2c — Extract LucaProt protein FAA

**Input:**
```
result/03_RDRP_identification/3_lucaprot/{sample}/RdRPs_only_using_threshold0.900000.csv
```

**Output:**
```
result/03_RDRP_identification/3_lucaprot/{sample}/{sample}_lucaprot_proteins.faa
```

ICTV equivalents follow same structure under `/scratch/xddeng/yangtze/RNA/my_rna/result/03_RDRP_identification/3_lucaprot/ICTV/`.

---

## Step 3 — Concatenate all candidate proteins  ← NEW

Combines all per-sample RdRpCATCH trimmed AA + LucaProt proteins into a single FAA
for the motif search in Step 4.

**Input:**
```
result/03_RDRP_identification/1_RdRpCATCH/{sample}/{sample}_scaffolds_rename_1000_rdrpcatch_fasta/
    {sample}_scaffolds_rename_1000_trimmed_aminoacid_sequences.fasta   (all samples, concatenated)

result/03_RDRP_identification/3_lucaprot/{sample}/{sample}_lucaprot_proteins.faa
    (all samples, concatenated)
```

**Output:**
```
result/03_RDRP_identification/9_motif_search/all_candidates.faa
```

This file is the single input FAA for all Step 4 motif search tools.

---

## Step 4 — RdRp motif A–D search (RVMT Sequence_Library)  ← NEW

Searches `all_candidates.faa` against RVMT motif profiles (mot.1–mot.4) using
five complementary tools to identify RdRp motif A, B, C, D.

Based on: Neri et al. (RVMT), RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library/

### Step 4a — Build profiles (once per motif)

**Input:**
```
database/RVMT/RVMT_Zenodo_V4/RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library/mot.{1,2,3,4}/
    *_woCon.afa   ← pre-aligned MSA cluster files
```

**Output (profile databases, per motif):**
```
result/03_RDRP_identification/9_motif_search/profiles/mot.{motif}/
    msaFiles/
        *.afa                    ← copied original MSAs
        *.Cons.msa.afa           ← consensus-added MSAs (hhconsensus -M 50 -cov 50)
        .consensus_done          ← sentinel file
    HMMfiles/
        *.hmm                    ← per-cluster HMMER3 profiles
        profiles.hmm             ← concatenated HMMER3 database (hmmbuild)
    HHMfiles/
        *.hhm                    ← per-cluster HH-suite profiles (hhmake)
        db / db.index            ← ffindex database (ffindex_build)
    MMseqs2_profiles/
        MM                       ← MMseqs2 profile database (mmseqs convertprofiledb)
    singletons.faa               ← MSA clusters with exactly 1 sequence (for Diamond)
    logs/
```

### Step 4b — Motif search (5 tools × 4 motifs = 20 output TSVs)

**Input for all search rules:**
```
result/03_RDRP_identification/9_motif_search/all_candidates.faa   ← your sequences
result/03_RDRP_identification/9_motif_search/profiles/mot.{motif}/ ← profile databases
```

**Output:**
```
result/03_RDRP_identification/9_motif_search/search/
    hhalign/mot.{1,2,3,4}.tsv      ← HHsearch profile-profile (highest sensitivity)
    hmmsearch/mot.{1,2,3,4}.tsv    ← HMMER3 hmmsearch (domtblout, # lines stripped)
    psiblast/mot.{1,2,3,4}.tsv     ← PSI-BLAST (motif MSA vs your seqs)
    mmseqs/mot.{1,2,3,4}.tsv       ← MMseqs2 profile search
    diamond/mot.{1,2,3,4}.tsv      ← Diamond BLASTp (singletons vs your seqs)
```

#### Output column formats

**hhalign** (13 cols, HHsearch blasttab):
```
subject_name  pCoverage  ali_len  pL  mismatch  gapOpen
q1  q2   ← start/end on YOUR sequence (AA)
p1  p2   ← start/end on motif profile
Probab  evalue  score
```

**hmmsearch** (23 cols, HMMER domtblout):
```
target_name(your_seq)  target_acc  qL
query_name(motif)  query_acc  pL
E-value  score  bias  #  of  c-Evalue  i-Evalue  score2  bias2
p1  p2   ← start/end on motif profile (cols 16-17)
q1  q2   ← start/end on YOUR sequence (cols 18-19)
env_from  env_to  acc  description
```

**psiblast** (12 cols, NOTE: query=motif MSA, subject=your seq):
```
sseqid(your_seq)  pident
sstart  send   ← start/end on YOUR sequence  (= q1, q2)
qstart  qend   ← start/end on motif MSA      (= p1, p2)
slen(your_seq_len)  qlen(motif_len)  ali_len  evalue  bitscore  profile_name
```

**mmseqs** (19 cols):
```
query(your_seq)  target(motif)  evalue  gapopen  pident  nident
qstart  qend  qlen   ← start/end/len on YOUR sequence
tstart  tend  tlen   ← start/end/len on motif profile
alnlen  raw  bits  qframe  mismatch  qcov  tcov
```

**diamond** (14 cols, NOTE: query=singleton, subject=your seq):
```
sseqid(your_seq)  qseqid(mot.{motif}.singleton_id)  evalue  gapopen  pident  ali_len
qstart  qend   ← start/end on singleton    (= p1, p2)
sstart  send   ← start/end on YOUR sequence (= q1, q2)
qlen(singleton_len)  slen(your_seq_len)  mismatch  bitscore
```

---

## Notes

- Steps 1–2 produce per-sample files; Step 3 merges them all for Step 4.
- Step 4 profiles (Step 4a) are built once and reused for all search tools.
- The `q1/q2` coordinates in the output TSVs are positions on YOUR sequence
  (all_candidates.faa). Use these for motif localisation and de-permutation.
- `psiblast` and `diamond` have reversed query/subject — see column notes above.
- Downstream parsing (R script `parse_motif_hits.R`) applies stricter E-value
  and coverage filters than the loose cutoffs used during search.
- ICTV sequences are NOT included in Step 3/4 — motif search is Contig-only.
  Run separately if needed.
