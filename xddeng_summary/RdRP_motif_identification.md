# RdRP Motif A–D Identification in RVMT

## Overview

Motifs A, B, C, and D are the four conserved palm-domain motifs of RNA-dependent RNA polymerases (RdRPs). In RVMT they are identified entirely through **profile-based similarity searches** — not by scanning for fixed amino acid patterns. Each motif is assigned an internal number (`mot.1` = A, `mot.2` = B, `mot.3` = C, `mot.4` = D).

The pipeline has two phases:
1. **Build profile databases** from known motif sequences (`Profiler_motifs.sh`)
2. **Search** candidate RdRP sequences against those databases (five parallel runner scripts)

---

## Phase 1: Building Motif Profile Databases

**Script:** `Domains_Annotation/Profiler_motifs.sh`

Run once per motif (e.g., `-P Vfin.mot.1` for motif A). Each run takes a FASTA file of known sequences covering that motif region and produces HMM, HH-suite, and MMseqs2 profile databases.

### Pipeline inside `Profiler_motifs.sh`

```
Input .faa (known motif-region sequences)
        │
        ▼
[1] Pre-cluster with MMseqs2 easy-linclust
    --min-seq-id $min_prec_id (default 0.95)
    --cov-mode 1, -c $min_prec_cov (default 0.95)
        │
        ▼
[2] All-vs-all Diamond BLASTp
    diamond blastp -q input --db input
    --query-cover 0.45 --subject-cover 0.45 -e 0.0001
        │
        ▼
[3] MCL clustering
    mcxload  →  mcl -I $MCL_inflation (default 1.4)
        │
        ▼
[4] Split input into per-cluster .faa files
    (seqkit grep or awk)
        │
        ▼
[5] Align each cluster
    MUSCLE5 or MAFFT (parallel)
        │
        ▼
[6] Add consensus sequence
    hhconsensus -M 50 -cov 50
        │
        ▼
[7] Build profile databases (parallel)
    HMMER  →  hmmbuild  →  HMMdb (.hmm)
    HH-suite  →  hhmake + ffindex  →  HHM db
    MMseqs2  →  mmseqs createdb + convertmsa + msa2profile  →  MMseqs2 profile db
```

**Key parameters (defaults in script):**

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `-d` min_prec_id | 0.95 | Min identity for pre-clustering |
| `-c` min_prec_cov | 0.95 | Min coverage for pre-clustering |
| `-f` MCL_inflation | 1.4 | MCL granularity (higher = finer clusters) |
| `-N` min_nseq | 10 | Min sequences per cluster to build a profile |
| `-L` Use_singlt | True | Include singleton clusters as profiles |

**Output directory structure:**
```
output/
├── msaFiles/          # per-cluster alignments (.msa.faa, .Cons.msa.afa)
├── singletons/        # single-sequence clusters
├── HMMfiles/          # HMMER3 profiles (.hmm)
├── HHMfiles/          # HH-suite profiles (.hhm)
├── MMseqs2_profiles/  # MMseqs2 profile database
└── singletons.faa     # merged singleton sequences
```

---

## Phase 2: Searching Candidate RdRPs Against Motif Profiles

Five runner scripts, each accepting `Motif_type` (1–4) as argument 4/5, search candidate sequences against the motif database and write results tagged with the motif number.

### 2a. HHsearch (profile–profile) — `Motif_hhseach_runner.sh`

```bash
# For each input profile alignment (Multi mode):
hhsearch -M 50 -e 0.5 \
    -i $input/{profile.afa} \
    -d $motifDB_path \
    -o /dev/null -cpu 1 -hide_cons -Z 100000 \
    -blasttab HHMsearch_vs_mot${Motif_type}/${seq}_hhsearch_rawout.tsv

# Concatenate all per-sequence results:
cat HHMsearch_*/* > ${out_pr}_HHMsearch_vs_${ini_name}_vs_mot${Motif_type}_${out_su}_rawout.tsv

# Clean up consensus suffix from profile names:
sed -i "s|.Cons.msa||g" output.tsv
```

**E-value cutoff:** 0.5 (loose; downstream R parsing applies stricter filters)

### 2b. HMMER hmmsearch — `Motif_hmmsearch_runner.sh`

```bash
hmmsearch --noali --cpu $THREADS \
    -E 0.5 --incE 0.5 \
    --domtblout ${out_pr}_HMMsearch_vs_${ini_name}_vs_mot${Motif_type}_${out_su}.tsv \
    $motifDB_path $input

# Strip comment lines, then reformat to clean TSV (23 fields):
sed -i '/^#/d' output.tsv
awk '{print $1"\t"$2...$23}' output.tsv > ../output.tsv
```

### 2c. PSI-BLAST — `Motif_psiblast_runner.sh`

```bash
# Iterates over every consensus MSA in the motif DB directory:
for cluster in "$motifDB_path"/*Cons.*.faa; do
    psiblast -word_size 2 -evalue 0.5 \
        -max_target_seqs 100000 -threshold 9 \
        -dbsize 20000000 -ignore_msa_master \
        -in_msa $cluster \
        -qcov_hsp_perc 1 \
        -num_threads $THREADS \
        -db $blastp_DB \
        -outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore" \
        -out ${cluster_name}_pisblast_${ini_name}.tsv
done

# Concatenate, appending profile name as final column:
ls *_pisblast_*.tsv | xargs -I% sed 's/$/\t%/' % > ${out_pr}_psiblast_vs_mot${Motif_type}.tsv

# Strip suffixes to leave clean profile cluster name:
sed -i "s|.Cons${tbf}||g" output.tsv
```

### 2d. MMseqs2 profile search — `Motif_MMseqs_runner.sh`

```bash
mmseqs createdb $input ${ini_name}_MMseqs2DB/MMseqs2DB

mmseqs search $querydb $motifDB_path $resultsdb ./tmp/ \
    -k 6 -s 7.5 -e 0.5 \
    --threads $THREADS --split-memory-limit 100G --max-seqs 6000

mmseqs convertalis $querydb $motifDB_path $resultsdb output.tsv \
    --format-output "query,target,evalue,gapopen,pident,nident,\
                     qstart,qend,qlen,tstart,tend,tlen,\
                     alnlen,raw,bits,qframe,mismatch,qcov,tcov"
```

### 2e. Diamond BLASTp — `Motif_DiamondP_runner.sh`

```bash
diamond makedb --in $input -d ${ini_name}_Diamond_DB/${ini_name}.dmnd

# Query = motif DB singletons; subject = candidate RdRP sequences
diamond blastp \
    -q $motifDB_path \
    --db $Diamond_DB \
    -b3.0 -k 0 \
    --query-cover 50 --id 60 --more-sensitive \
    --evalue 0.05 \
    --outfmt 6 sseqid qseqid evalue gapopen pident length \
                qstart qend sstart send qlen slen mismatch bitscore \
    -o ${out_pr}_DiamondP_vs_${ini_name}_vs_mot${Motif_type}_${out_su}.tsv

# Prepend motif type label (e.g. "mot.1.") to the query (profile) column:
awk -F'\t' -vOFS='\t' -v Motif_type="$Motif_type" \
    '{ $2 = "mot."Motif_type"."$2 }1' output.tsv > tmpout && mv tmpout output.tsv
```

---

## Phase 3: Parsing Hit Positions in R

**Function:** `GenericHitsParser()` in `Utils/basicf.r` (line 919)

Standardizes all five search output formats into a unified data frame with these key columns:

| Column | Meaning |
|--------|---------|
| `q1`, `q2` | Start/end of the motif match on the **candidate RdRP sequence** (AA coordinates) |
| `p1`, `p2` | Start/end of the match on the **motif profile** |
| `qL` | Full length of the candidate sequence |
| `pL` | Length of the motif profile |
| `score`, `evalue` | Hit quality |
| `subject_name` | Which motif profile cluster was matched (includes `mot.1` – `mot.4` label) |

Column name sets used by the parser:

```r
# Defined in Utils/basicf.r lines 908–913
MMseq2_outfmt6_cols  <- c("subject_name","evalue","gapopen","pident","nident",
                           "q1","q2","qL","p1","p2","pL","ali_len","raw","score",
                           "frame","mismatch","qcov","tcov")

hmmsearch_cols       <- c("r1","qL","subject_name","r2","pL","evalue","score",
                           "bias","#","of","c-Evalue","i-Evalue","score2","bias",
                           "p1","p2","q1","q2","env_from","env_to","acc","r3","r4")

psiblast_cols        <- c("pident","q1","q2","p1","p2","qL","pL",
                           "ali_len","evalue","score","subject_name")

hhsearch_cols        <- c("subject_name","pCoverage","ali_len","pL","mismatch",
                           "gapOpen","q1","q2","p1","p2","Probab","evalue","score")

DaimondP_cols        <- c("subject_name","evalue","gapopen","pident","ali_len",
                           "p1","p2","q1","q2","pL","qL","mismatch","score")
```

**Post-parse filters applied in `Annotations.r`:**

```r
# Example thresholds from the RNAVirDB search block (Annotations.r lines 78–81):
data.table(HRVA)[evalue  <= 1e-5]
data.table(HRVA)[ali_Qcov_len >= 25]
data.table(HRVA)[score   >= 10]
filter(HRVA, !(PolyPort_Problematic))
```

**Sequence extraction from hit coordinates** (`Trim2Core2()`, `Utils/basicf.r` line 875):

```r
Trim2Core2 <- function(hits_df, faa, Expand_projectionX = 0, ...) {
    # Cuts the motif region out of each full-length RdRP sequence
    # using q1/q2 from the search output, optionally expanding by X residues
    out_faa <- BStringSet(hits_df$seq, start = hits_df$q1, end = hits_df$q2)
    # With expansion:
    # start = max(1,  q1 - Expand_projectionX * (p1 - 1))
    # end   = min(qL, q2 + Expand_projectionX * (pL - p2))
}
```

---

## Motif Order and Permutation Detection

After parsing, each candidate RdRP has `q1` coordinates recorded for all four motifs. The **canonical order** is motif A → B → C → D (ascending `q1`). Sequences where motif C (`mot.3`) has a lower `q1` than motifs A and B are flagged as **permuted** (C-A-B-D configuration).

The permutation status is stored in the `RBS` field of the master metadata table (`IDFT`) and in `Perm_Clades.txt`. These feed into the tree annotation scripts:

```r
# Phylogeny/Tree_Plots.R lines 182–190
PermutSets <- fread("Perm_Clades.txt", sep = "\t")
PermutSets$leaflist <- apply(PermutSets, 1, function(x)
    as.char(unlist(str_split_fixed(x["V2"], ",", n = Inf))))
PermutSets$node       <- unlist(apply(PermutSets, 1, function(x)
    get_mrca_of_set(WorkTree, unlist(x["leaflist"]))))
PermutSets$node.label <- NID2NLabel(WorkTree, PermutSets$node)

# Mark permuted clades on tree:
info2$shapo[wh(info2$node.label %in% PermutSets$node.label)] <- "Permuted_RdRP"
```

```r
# Phylogeny/Tree_Plots.R lines 335–339
p12$data$Perm[wh(p12$data$node %in% PermutSets$node)] <- TRUE
p12$data$shapo[wh(p12$data$Perm)] <- "Permuted RdRP"
```

The actual per-sequence comparison of `q1(mot.3) < q1(mot.1)` to call permutation was done outside this repository; only the resulting `Perm_Clades.txt` and the `RBS` column are retained here.

---

## Summary of Search Parameters

| Script | Tool | E-value | Query cover | Identity |
|--------|------|---------|------------|---------|
| `Motif_hhseach_runner.sh` | HHsearch | 0.5 | — (profile–profile) | — |
| `Motif_hmmsearch_runner.sh` | HMMER hmmsearch | 0.5 | — | — |
| `Motif_psiblast_runner.sh` | PSI-BLAST | 0.5 | 1% | — |
| `Motif_MMseqs_runner.sh` | MMseqs2 | 0.5 | — | — |
| `Motif_DiamondP_runner.sh` | Diamond BLASTp | 0.05 | 50% | ≥60% |

All searches use loose E-value cutoffs at the search stage; strict filtering (`evalue ≤ 1e-5`, `score ≥ 10`, `ali_len ≥ 25`) is applied during R parsing in `Annotations.r`.
