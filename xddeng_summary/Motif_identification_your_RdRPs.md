# Identifying RdRP Motifs A–D in Your Own Sequences

## What You Already Have

| Path | Content |
|------|---------|
| `RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library/mot.1–4/` | Pre-built motif cluster MSAs (`.afa` files, ~1173 clusters total) |
| Your candidate RdRP sequences | Full-length or domain-trimmed protein sequences in `.faa` format |

You **do not** need to re-cluster or re-align anything. The `Sequence_Library` MSAs are the direct input to profile building, and the pipeline skips directly to that step using `-k True` in `Profiler_motifs.sh`.

---

## Overview of Steps

```
Your candidate RdRPs (.faa)
        │
        ▼
[Step 1] Build profile databases from Sequence_Library MSAs
         (Profiler_motifs.sh -k True, run once per motif)
         → HMM / HH-suite / MMseqs2 profile DBs  (mot.1 – mot.4)
        │
        ▼
[Step 2] Search your RdRPs against each motif profile DB
         (run each runner script 4×, once per motif)
         → Raw hit tables (.tsv)
        │
        ▼
[Step 3] Parse and filter hits in R
         (GenericHitsParser + thresholds)
         → Unified hit table with q1, q2, motif_type per sequence
        │
        ▼
[Step 4] Determine motif order per sequence
         (compare q1 across mot.1–4)
         → Canonical (A-B-C-D) or Permuted (C-A-B-D) flag
```

---

## Step 1: Build Profile Databases from the Sequence Library

**Script:** `Domains_Annotation/Profiler_motifs.sh`

Run **four times** — once per motif. The `-k True` flag skips clustering and alignment; it treats the input directory as already-aligned MSAs.

### Input

```
RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library/mot.1/
├── mot.1.22.1_woCon.afa      # aligned FASTA, no consensus sequence
├── mot.1.22.2_woCon.afa
├── mot.1.25.0.1_woCon.afa
└── ...                        # 296 files total for mot.1
```

Each `.afa` file is a multi-sequence alignment in FASTA format. Sequence IDs follow the format `SeqID/start-end`:
```
>1134955769/1-24
G-AKFYSVDLSAASDRLSQPLSLSV
>1134599847/1-23
P-WVALSCDMRSATDNFPHYLVE-A
```

### Command (run 4×, change `mot.1` → `mot.2`, `mot.3`, `mot.4`)

```bash
bash Profiler_motifs.sh \
    -t 16 \
    -M 50000 \
    -o /output/motif_profiles/mot.1/ \
    -i /RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library/mot.1/ \
    -P Vfin.mot.1 \
    -k True \
    -a True \
    -L True \
    -N 2 \
    -r False
```

Key flags:
| Flag | Value | Meaning |
|------|-------|---------|
| `-k True` | True | **Skip** clustering/alignment; use existing MSAs |
| `-i` | path to `mot.X/` directory | Input directory of `.afa` files |
| `-P` | `Vfin.mot.1` | Output cluster name prefix |
| `-a True` | True | Add consensus sequence to each alignment |
| `-L True` | True | Include singleton clusters |

### Output

```
/output/motif_profiles/mot.1/
├── msaFiles/
│   ├── Vfin.mot.1.1.msa.faa          # per-cluster alignment
│   ├── Vfin.mot.1.1.Cons.msa.afa     # alignment + consensus on top
│   └── ...
├── HMMfiles/
│   └── Vfin.mot.1_profiles.hmm       # HMMER3 profile database
├── HHMfiles/
│   ├── db                             # HH-suite profile database
│   └── db.index
├── MMseqs2_profiles/
│   └── MMseqs2DB                      # MMseqs2 profile database
└── singletons.faa                     # merged singleton sequences
```

---

## Step 2: Search Your RdRPs Against Each Motif Profile Database

Run each search tool 4× (once per motif). All scripts share the same argument structure:
`THREADS  output_dir  input.faa  Motif_type  motifDB_path  out_prefix  out_suffix  rm_tmp`

### Your input sequences

```
your_RdRPs.faa
>ND_432610..3..1:300
QLLSAWKAHPTSTMTVNLDEAVQKLY...
>ND_432609..2..220:577
ELNTLVNRGYGAVNWQTEREHRLNPD...
...
```

Standard protein FASTA. Sequence IDs can be anything; the parsers split on `.` to extract contig ID and frame.

---

### 2a. HHsearch (profile–profile, highest sensitivity)

**Script:** `Discovery_pipeline/RdRP_searchs/Motif_hhseach_runner.sh`

```bash
bash Motif_hhseach_runner.sh \
    16 \                                          # $1 THREADS
    /output/hhsearch_mot1/ \                      # $2 output dir
    single \                                      # $3 input type
    your_RdRPs.faa \                              # $4 input
    1 \                                           # $5 Motif_type (1–4)
    /output/motif_profiles/mot.1/HHMfiles/db \   # $6 motifDB path
    MyRun \                                       # $7 out prefix
    v1 \                                          # $8 out suffix
    False                                         # $9 remove tmps
```

**What it does internally:**
```bash
# Splits your_RdRPs.faa into one file per sequence
splitfasta.pl your_RdRPs.faa -ext .faa

# Runs HHsearch in parallel, one job per sequence
parallel hhsearch -M 50 -e 0.5 \
    -i {seq}.faa \
    -d /output/motif_profiles/mot.1/HHMfiles/db \
    -o /dev/null -cpu 1 -hide_cons -Z 100000 \
    -blasttab HHMsearch_vs_mot1/{seq}_hhsearch_rawout.tsv

# Concatenates all results
cat HHMsearch_vs_mot1/* > MyRun_HHMsearch_vs_mot1_v1_rawout.tsv

# Cleans up consensus suffix from profile names
sed -i "s|.Cons.msa||g" MyRun_HHMsearch_vs_mot1_v1_rawout.tsv
```

**Output:** `MyRun_HHMsearch_vs_mot1_v1_rawout.tsv`

```
# Tab-separated, columns (hhsearch blasttab format):
subject_name    pCoverage   ali_len   pL   mismatch   gapOpen
q1   q2   p1   p2   Probab   evalue   score
```

Example row:
```
Vfin.mot.1.22   0.96   24   26   1   0   194   217   1   25   99.2   3.1e-09   48.3
```

---

### 2b. HMMER hmmsearch

**Script:** `Discovery_pipeline/RdRP_searchs/Motif_hmmsearch_runner.sh`

```bash
bash Motif_hmmsearch_runner.sh \
    16 \
    /output/hmmsearch_mot1/ \
    your_RdRPs.faa \
    1 \
    /output/motif_profiles/mot.1/HMMfiles/Vfin.mot.1_profiles.hmm \
    MyRun \
    v1 \
    False
```

**What it does internally:**
```bash
hmmsearch --noali --cpu 16 -E 0.5 --incE 0.5 \
    --domtblout MyRun_HMMsearch_vs_mot1_v1.tsv \
    /output/motif_profiles/mot.1/HMMfiles/Vfin.mot.1_profiles.hmm \
    your_RdRPs.faa

# Remove comment lines
sed -i '/^#/d' MyRun_HMMsearch_vs_mot1_v1.tsv

# Reformat to 23-column TSV
awk '{print $1"\t"$2"\t"...$23}' MyRun_HMMsearch_vs_mot1_v1.tsv > ../MyRun_HMMsearch_vs_mot1_v1.tsv
```

**Output:** `MyRun_HMMsearch_vs_mot1_v1.tsv`

```
# 23 columns (hmmer domtblout format, comments stripped):
target_name   target_accession   qL   query_name   query_accession
pL   E-value   score   bias   #   of   c-Evalue   i-Evalue
bias2   score2   p1   p2   q1   q2   env_from   env_to   acc   description
```

---

### 2c. PSI-BLAST

**Script:** `Discovery_pipeline/RdRP_searchs/Motif_psiblast_runner.sh`

```bash
bash Motif_psiblast_runner.sh \
    16 \
    /output/psiblast_mot1/ \
    your_RdRPs.faa \
    1 \
    /output/motif_profiles/mot.1/msaFiles/ \
    MyRun \
    v1 \
    False
```

**What it does internally:**
```bash
# Build a BLAST database from your RdRPs
makeblastdb -in your_RdRPs.faa -dbtype prot -out your_RdRPs_DB

# Search each motif cluster consensus MSA against your RdRPs
for cluster in /output/motif_profiles/mot.1/msaFiles/*Cons.*.faa; do
    psiblast -word_size 2 -evalue 0.5 \
        -max_target_seqs 100000 -threshold 9 \
        -dbsize 20000000 -ignore_msa_master \
        -in_msa $cluster \
        -qcov_hsp_perc 1 \
        -num_threads 16 \
        -db your_RdRPs_DB \
        -outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore" \
        -out ${cluster_name}_psiblast.tsv
done

# Concatenate, appending cluster name as last column
ls *_pisblast_*.tsv | xargs -I% sed 's/$/\t%/' % > MyRun_psiblast_vs_mot1_v1.tsv

# Strip suffixes to leave clean cluster name
sed -i "s|.Cons_pisblast_your_RdRPs.tsv||g" MyRun_psiblast_vs_mot1_v1.tsv
```

**Output:** `MyRun_psiblast_vs_mot1_v1.tsv`

```
# 12 columns:
sseqid(=your RdRP ID)   pident   sstart   send   qstart   qend
slen   qlen   length   evalue   bitscore   profile_cluster_name
```

---

### 2d. MMseqs2 profile search

**Script:** `Discovery_pipeline/RdRP_searchs/Motif_MMseqs_runner.sh`

```bash
bash Motif_MMseqs_runner.sh \
    16 \
    /output/mmseqs_mot1/ \
    your_RdRPs.faa \
    1 \
    /output/motif_profiles/mot.1/MMseqs2_profiles/MMseqs2DB \
    MyRun \
    v1 \
    False
```

**What it does internally:**
```bash
mmseqs createdb your_RdRPs.faa your_RdRPs_MMseqs2DB/MMseqs2DB

mmseqs search your_RdRPs_MMseqs2DB/MMseqs2DB \
    /output/motif_profiles/mot.1/MMseqs2_profiles/MMseqs2DB \
    resultsdb ./tmp/ \
    -k 6 -s 7.5 -e 0.5 \
    --threads 16 --split-memory-limit 100G --max-seqs 6000

mmseqs convertalis \
    your_RdRPs_MMseqs2DB/MMseqs2DB \
    /output/motif_profiles/mot.1/MMseqs2_profiles/MMseqs2DB \
    resultsdb \
    MyRun_MMseqs2_vs_mot1_v1.tsv \
    --format-output "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,qframe,mismatch,qcov,tcov"
```

**Output:** `MyRun_MMseqs2_vs_mot1_v1.tsv`

```
# 19 columns:
query   target   evalue   gapopen   pident   nident
qstart(=q1)   qend(=q2)   qlen(=qL)
tstart(=p1)   tend(=p2)   tlen(=pL)
alnlen   raw   bits   qframe   mismatch   qcov   tcov
```

---

### 2e. Diamond BLASTp

**Script:** `Discovery_pipeline/RdRP_searchs/Motif_DiamondP_runner.sh`

```bash
bash Motif_DiamondP_runner.sh \
    16 \
    /output/diamond_mot1/ \
    your_RdRPs.faa \
    1 \
    /output/motif_profiles/mot.1/singletons.faa \
    MyRun \
    v1 \
    False
```

**What it does internally:**
```bash
# Build Diamond DB from your RdRPs
diamond makedb --in your_RdRPs.faa -d your_RdRPs_Diamond_DB/your_RdRPs.dmnd

# Query = motif singleton sequences; subject = your RdRPs
diamond blastp \
    -q /output/motif_profiles/mot.1/singletons.faa \
    --db your_RdRPs_Diamond_DB/your_RdRPs.dmnd \
    -b3.0 -k 0 \
    --query-cover 50 --id 60 --more-sensitive \
    --evalue 0.05 \
    --outfmt 6 sseqid qseqid evalue gapopen pident length \
                qstart qend sstart send qlen slen mismatch bitscore \
    -o MyRun_DiamondP_vs_mot1_v1.tsv

# Prepend "mot.1." to the query (motif profile) column
awk -F'\t' -vOFS='\t' -v Motif_type="1" \
    '{ $2 = "mot."Motif_type"."$2 }1' \
    MyRun_DiamondP_vs_mot1_v1.tsv > tmpout && mv tmpout MyRun_DiamondP_vs_mot1_v1.tsv
```

**Output:** `MyRun_DiamondP_vs_mot1_v1.tsv`

```
# 14 columns:
sseqid(=your RdRP ID)   qseqid(=mot.1.cluster_name)
evalue   gapopen   pident   length
qstart(=p1)   qend(=p2)   sstart(=q1)   send(=q2)
qlen(=pL)   slen(=qL)   mismatch   bitscore
```

---

## Step 3: Parse and Filter in R

All 20 raw TSV files (5 tools × 4 motifs) are loaded with `GenericHitsParser()` from `Utils/basicf.r`,
merged into one table, then filtered. A standalone version of this script is at
[parse_motif_hits.R](parse_motif_hits.R).

### Column name vectors (defined in `Utils/basicf.r` lines 908–916)

These tell `GenericHitsParser()` how to interpret each tool's raw TSV:

```r
# HHsearch blasttab (13 cols)
hhsearch_cols <- c("subject_name","pCoverage","ali_len","pL","mismatch",
                   "gapOpen","q1","q2","p1","p2","Probab","evalue","score")

# HMMER domtblout (23 cols, comments stripped)
hmmsearh_cols <- c("r1","qL","subject_name","r2","pL","evalue","score",
                   "bias","#","of","c-Evalue","i-Evalue","score2","bias",
                   "p1","p2","q1","q2","env_from","env_to","acc","r3","r4")

# PSI-BLAST outfmt6 + profile name appended (12 cols)
psiblast_cols <- c("pident","q1","q2","p1","p2","qL","pL",
                   "ali_len","evalue","score","subject_name")

# MMseqs2 convertalis (19 cols, nident/frame/mismatch dropped by parser)
MMseq2_outfmt6_cols <- c("subject_name","evalue","gapopen","pident","nident",
                          "q1","q2","qL","p1","p2","pL","ali_len","raw",
                          "score","frame","mismatch","qcov","tcov")

# Diamond BLASTp (14 cols; q1/q2 are sstart/send because query/subject reversed)
DaimondP_cols <- c("subject_name","evalue","gapopen","pident","ali_len",
                   "p1","p2","q1","q2","pL","qL","mismatch","score")
```

### `GenericHitsParser()` — key arguments

```r
GenericHitsParser(
    inpt          = "path/to/file.tsv",  # raw search output file
    search_tool   = "hhsearch",          # one of: hhsearch, hmmsearch, psiblast, mmseqs, diamondp
    input_was     = "RdRp_id",           # column name to assign to the first column (your seq ID)
    Query2Profile = TRUE,                # TRUE = your seq is the query (all tools except psiblast/diamond)
    breakhdrs     = FALSE,               # FALSE = keep sequence IDs as-is
    reducecols    = TRUE,                # keep only unified columns
    calc_pcoverage = TRUE                # compute pCoverage = (p2-p1+1)/pL
)
```

`GenericHitsParser()` internally:

1. Reads the TSV with `fread()`, assigning the correct column names for the tool
2. For `mmseqs`: drops columns 6 (nident), 4 (gapopen), 16 (frame)
3. For `hmmsearch`: drops redundant columns and computes `ali_len = q2 - q1`
4. Renames `subject_name` → `profile`
5. Computes `pCoverage = (p2 - p1 + 1) / pL` if `calc_pcoverage = TRUE`
6. Returns a data frame with unified columns: `RdRp_id, profile, qL, pL, p1, p2, q1, q2, score, evalue, ali_len, pCoverage`

### Full parse and merge loop

```r
library(data.table)
library(Biostrings)
source("Utils/basicf.r")   # provides GenericHitsParser(), Trim2Core2(), CalcPcoverage()

OUTDIR <- "/path/to/output"   # must match config.yaml outdir
INPUT_FAA <- "your_RdRPs.faa"

motifs <- 1:4
tools  <- c("hhalign", "hmmsearch", "psiblast", "mmseqs", "diamond")

# Map tool folder names to GenericHitsParser search_tool argument
tool_map <- c(
    hhalign   = "hhsearch",
    hmmsearch = "hmmsearch",
    psiblast  = "psiblast",
    mmseqs    = "mmseqs",
    diamond   = "diamondp"
)

all_hits <- rbindlist(lapply(motifs, function(m) {
    rbindlist(lapply(tools, function(tool) {
        f <- file.path(OUTDIR, "search", tool, paste0("mot.", m, ".tsv"))
        if (!file.exists(f) || file.size(f) == 0) return(NULL)

        h <- GenericHitsParser(
            inpt          = f,
            search_tool   = tool_map[tool],
            input_was     = "RdRp_id",
            Query2Profile = TRUE,
            breakhdrs     = FALSE,
            reducecols    = TRUE,
            calc_pcoverage = TRUE
        )
        h$motif_type <- m
        h$search_tool <- tool
        h
    }), fill = TRUE)
}), fill = TRUE)
```

### Filter thresholds (from `Annotations.r` lines 907–909)

```r
# Convert types (fread may read as character)
all_hits$evalue <- as.numeric(all_hits$evalue)
all_hits$score  <- as.numeric(all_hits$score)
all_hits$q1     <- as.integer(all_hits$q1)
all_hits$q2     <- as.integer(all_hits$q2)
all_hits$ali_len <- as.integer(all_hits$q2 - all_hits$q1)

# Apply Neri's primary filter thresholds
hits_filt <- all_hits[evalue <  1e-5]
hits_filt <- hits_filt[score  >  10  ]
hits_filt <- hits_filt[ali_len > 10  ]

# Cull to best hit per sequence per motif (across all tools)
setorder(hits_filt, RdRp_id, motif_type, -score, evalue)
hits_best <- hits_filt[, .SD[1], by = .(RdRp_id, motif_type)]
```

### Unified output table (mirrors `ite012_unfiltered.tsv`)

| Column | Example | Meaning |
| ------ | ------- | ------- |
| `RdRp_id` | `seq_001` | Your sequence ID |
| `motif_type` | `1` | Motif number (1=A, 2=B, 3=C, 4=D) |
| `q1` | `193` | Start of motif on your sequence (AA) |
| `q2` | `217` | End of motif on your sequence (AA) |
| `qL` | `427` | Your sequence length |
| `p1` | `2` | Start on motif profile |
| `p2` | `26` | End on motif profile |
| `pL` | `26` | Motif profile length |
| `score` | `37.2` | Bit score |
| `evalue` | `4.4e-6` | E-value |
| `pCoverage` | `0.96` | Profile coverage = (p2−p1+1)/pL |
| `profile` | `Vfin.mot.1.22` | Matched motif cluster |
| `ali_len` | `24` | Alignment length on your sequence |
| `search_tool` | `hmmsearch` | Which tool found this hit |

```r
# Write unified filtered table
fwrite(hits_best, "motif_hits_filtered.tsv", sep = "\t")
```

### Extract motif subsequences with `Trim2Core2()`

```r
# Trim2Core2() defined in Utils/basicf.r line 875
# Cuts q1:q2 from each full-length sequence
your_faa <- readAAStringSet(INPUT_FAA)

motif_seqs <- Trim2Core2(
    hits_df          = hits_best,
    faa              = your_faa,
    input_was        = "RdRp_id",   # column matching names(your_faa)
    subject_was      = "profile",   # used to build output sequence name
    Expand_projectionX = 0,         # 0 = exact q1:q2, no expansion
    Out_faa          = TRUE,
    out_df           = FALSE
)
# Returns: AAStringSet
# Sequence names: "SeqID.ProfileName"  e.g. "seq_001.Vfin.mot.1.22"
writeXStringSet(motif_seqs, "motif_sequences.faa")
```

---

## Step 4: Determine Motif Order (Canonical vs. Permuted)

For each sequence, compare `q1` positions across the four motifs to detect the
permuted RdRP configuration (C–A–B–D instead of A–B–C–D).

```r
library(tidyr)

# Pivot to wide format: one row per sequence, one column per motif
motif_pos <- dcast(
    hits_best,
    RdRp_id ~ motif_type,
    value.var = "q1",
    fun.aggregate = min   # take earliest start if multiple hits
)
setnames(motif_pos,
         old = c("1","2","3","4"),
         new = c("q1_A","q1_B","q1_C","q1_D"))

# How many motifs were found per sequence?
motif_pos$n_motifs_found <- rowSums(!is.na(motif_pos[, .(q1_A, q1_B, q1_C, q1_D)]))

# Canonical: A < B < C < D
motif_pos$canonical <- with(motif_pos,
    !is.na(q1_A) & !is.na(q1_B) & !is.na(q1_C) & !is.na(q1_D) &
    q1_A < q1_B & q1_B < q1_C & q1_C < q1_D
)

# Permuted: motif C is upstream of A  (C–A–B–D configuration)
motif_pos$permuted <- with(motif_pos,
    !is.na(q1_C) & !is.na(q1_A) & q1_C < q1_A
)

# Incomplete: fewer than 4 motifs found
motif_pos$incomplete <- motif_pos$n_motifs_found < 4

fwrite(motif_pos, "motif_order.tsv", sep = "\t")
```

### Output table

| RdRp_id | q1_A | q1_B | q1_C | q1_D | n_motifs_found | canonical | permuted | incomplete |
| ------- | ---- | ---- | ---- | ---- | -------------- | --------- | -------- | ---------- |
| seq_001 | 45   | 120  | 210  | 290  | 4              | TRUE      | FALSE    | FALSE      |
| seq_002 | 95   | 155  | 55   | 260  | 4              | FALSE     | TRUE     | FALSE      |
| seq_003 | NA   | 88   | 170  | 240  | 3              | FALSE     | FALSE    | TRUE       |

**Interpretation:**

- `canonical = TRUE` → sequence is ready for phylogenetic analysis as-is
- `permuted = TRUE` → motif C is upstream of A (C–A–B–D); sequence must be de-permuted before alignment
- `incomplete = TRUE` → fewer than 4 motifs detected; too fragmentary for tree inclusion

De-permutation means cutting the N-terminal C-motif region and reinserting it between B and D, then re-running this table to confirm canonical order.

---

## Complete File Flow Summary

```
Sequence_Library/mot.1–4/ (.afa MSAs)
    │
    │  Profiler_motifs.sh -k True  (×4)
    ▼
motif_profiles/mot.1–4/
    ├── HMMfiles/*.hmm
    ├── HHMfiles/db
    ├── MMseqs2_profiles/MMseqs2DB
    └── msaFiles/*Cons.*.faa  +  singletons.faa
    │
    │  Motif_*_runner.sh  (×5 tools × 4 motifs = 20 runs)
    ▼
Raw hit tables (mot1..mot4 × hhsearch/hmm/psiblast/mmseqs/diamond)
    *.tsv  —  columns vary by tool (q1, q2, p1, p2, score, evalue)
    │
    │  GenericHitsParser() + filters in R
    ▼
Unified hits table
    RdRp_id | motif_type | q1 | q2 | qL | score | evalue | profile | ...
    │
    │  Trim2Core2()
    ▼
Motif sequences (.faa)
    >SeqID.ProfileName
    GLATLDLRGASNSVFVEFVRSVIPP
    │
    │  pivot_wider + permuted flag
    ▼
Motif order table
    RdRp_id | q1_mot1 | q1_mot2 | q1_mot3 | q1_mot4 | permuted | canonical
```

## Search Parameter Reference

| Script | Tool | E-value | Query cover | Identity | Notes |
|--------|------|---------|-------------|----------|-------|
| `Motif_hhseach_runner.sh` | HHsearch | 0.5 | — | — | Profile–profile; most sensitive |
| `Motif_hmmsearch_runner.sh` | HMMER | 0.5 | — | — | HMM–sequence |
| `Motif_psiblast_runner.sh` | PSI-BLAST | 0.5 | 1% | — | MSA-guided BLAST |
| `Motif_MMseqs_runner.sh` | MMseqs2 | 0.5 | — | — | `-k 6 -s 7.5` (high sensitivity) |
| `Motif_DiamondP_runner.sh` | Diamond | 0.05 | 50% | ≥60% | Strictest at search stage |
| R parsing | — | ≤1e-5 | — | — | Applied after merging all tools |
