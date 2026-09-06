# Script Argument Reference — RdRP Motif Search Scripts

---

## 1. `Profiler_motifs.sh`
**Purpose:** Build profile databases from MSA files (run once per motif)
**Argument style:** Named flags (`-flag value`)

```bash
bash Profiler_motifs.sh \
    -t  THREADS          \  # number of CPU threads [default: 12]
    -M  MEMORY_MB        \  # memory in MB [default: 50000]
    -o  output_dir/      \  # output directory [default: pwd]
    -i  input/           \  # input: path to .faa file OR dir of .afa MSAs (if -k True)
    -P  Cls_Prefix       \  # output cluster name prefix, e.g. "Vfin.mot.1"
    -k  True/False       \  # skip clustering & alignment? True = use existing MSAs [default: False]
    -a  True/False       \  # add consensus sequence to each alignment? [default: True]
    -L  True/False       \  # include singleton clusters as profiles? [default: True]
    -N  min_nseq         \  # min sequences per cluster to build a profile [default: 10]
    -d  min_prec_id      \  # min identity for pre-clustering [default: 0.95]
    -c  min_prec_cov     \  # min coverage for pre-clustering [default: 0.95]
    -f  MCL_inflation    \  # MCL inflation value [default: 1.4]
    -p  True/False       \  # precluster with MMseqs2? [default: True]
    -x  True/False       \  # Diamond max-sensitivity mode? [default: False]
    -S  True/False       \  # use MAFFT instead of MUSCLE? [default: True]
    -H  True/False       \  # add secondary structure annotation? [default: False]
    -r  True/False       \  # remove temp files? [default: False]
    -s  "params_string"     # watermark string written to .env log file
```

**For your use case (existing MSAs, skip clustering):**
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

---

## 2. `Motif_hhalign_runner.sh`
**Purpose:** HHsearch profile–profile search of your RdRPs vs motif HH-suite DB
**Argument style:** Named flags (`-flag value`)

```bash
bash Motif_hhalign_runner.sh \
    -t  THREADS          \  # $THREADS  — number of CPU threads [default: 2]
    -o  output_dir/      \  # output directory [default: pwd]
    -y  input_type       \  # "single" (one .faa file) or "Multi" (dir of .afa alignments) [default: single]
    -i  input            \  # path to your_RdRPs.faa (single) or dir of .afa files (Multi)
    -M  Motif_type       \  # motif number/ID, e.g. 1, 2, 3, 4  [default: Mot]
    -d  motifDB_path     \  # path to HH-suite profile DB (the "db" file, no extension)
    -p  out_prefix       \  # prefix for output .tsv filename [default: prefff]
    -s  out_suffix       \  # suffix for output .tsv filename [default: sufff]
    -e  evalue           \  # E-value cutoff [default: 0.001]
    -r  True/False          # remove temp files? [default: False]
```

**Output file:**
```
{out_prefix}_HHMsearch_vs_{ini_name}_vs_mot{Motif_type}_{out_suffix}_rawout.tsv
```

**Example (run 4×, change `-M` and `-d` each time):**
```bash
bash Motif_hhalign_runner.sh \
    -t 16 \
    -o /output/hhalign_mot1/ \
    -y single \
    -i your_RdRPs.faa \
    -M 1 \
    -d /output/motif_profiles/mot.1/HHMfiles/db \
    -p MyRun \
    -s v1 \
    -e 0.5 \
    -r False
```

**Output columns (HHsearch blasttab format, 13 columns):**
```
col 1   subject_name     matched motif profile cluster name
col 2   pCoverage        profile coverage (p2-p1+1)/pL
col 3   ali_len          alignment length
col 4   pL               motif profile length
col 5   mismatch         number of mismatches
col 6   gapOpen          number of gap openings
col 7   q1               start on YOUR sequence (AA position)
col 8   q2               end on YOUR sequence (AA position)
col 9   p1               start on motif profile
col 10  p2               end on motif profile
col 11  Probab           HHsearch probability (0–100)
col 12  evalue           E-value
col 13  score            bit score
```

---

## 3. `Motif_hhseach_runner.sh`
**Purpose:** HHsearch using positional arguments (alternative to `hhalign_runner`)
**Argument style:** Positional (`$1 $2 ...`)

```bash
bash Motif_hhseach_runner.sh \
    $1  THREADS          \  # number of CPU threads
    $2  output_dir/      \  # output directory
    $3  input_type       \  # "single" or "Multi"
    $4  input            \  # path to your_RdRPs.faa or dir of .afa files
    $5  Motif_type       \  # motif number: 1, 2, 3, or 4
    $6  motifDB_path     \  # path to HH-suite profile DB
    $7  out_prefix       \  # prefix for output .tsv filename
    $8  out_suffix       \  # suffix for output .tsv filename
    $9  rm_tmp           \  # "True" or "False" — remove temp files
    $10 extension           # file extension: ".hhm" or leave blank for .afa
```

**Output file:**
```
{out_prefix}_HHMsearch_vs_{ini_name}_vs_mot{Motif_type}_{out_suffix}_rawout.tsv
```

**Example:**
```bash
bash Motif_hhseach_runner.sh \
    16 \
    /output/hhseach_mot1/ \
    single \
    your_RdRPs.faa \
    1 \
    /output/motif_profiles/mot.1/HHMfiles/db \
    MyRun \
    v1 \
    False \
    .afa
```

**Output columns:** Same 13-column HHsearch blasttab format as `hhalign_runner` above.

---

## 4. `Motif_hmmsearch_runner.sh`
**Purpose:** HMMER hmmsearch of your RdRPs vs motif HMM profile DB
**Argument style:** Positional (`$1 $2 ...`)

```bash
bash Motif_hmmsearch_runner.sh \
    $1  THREADS          \  # number of CPU threads
    $2  output_dir/      \  # output directory
    $3  input            \  # path to your_RdRPs.faa
    $4  Motif_type       \  # motif number: 1, 2, 3, or 4
    $5  motifDB_path     \  # path to .hmm profile database file
    $6  out_prefix       \  # prefix for output .tsv filename
    $7  out_suffix       \  # suffix for output .tsv filename
    $8  rm_tmp              # "True" or "False" — remove temp files
```

**Output file:**
```
{out_prefix}_HMMsearch_vs_{ini_name}_vs_mot{Motif_type}_{out_suffix}.tsv
```

**Example:**
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

**Output columns (HMMER domtblout format, 23 columns, comments stripped):**
```
col 1   target_name      your sequence ID
col 2   target_accession "-" (unused)
col 3   qL               length of your sequence (AA)
col 4   query_name       motif profile cluster name
col 5   query_accession  "-" (unused)
col 6   pL               motif profile length
col 7   E-value          full-sequence E-value
col 8   score            full-sequence bit score
col 9   bias             full-sequence bias
col 10  #                domain number (this domain)
col 11  of               total domains found
col 12  c-Evalue         conditional E-value (this domain)
col 13  i-Evalue         independent E-value (this domain)
col 14  score2           domain bit score
col 15  bias2            domain bias
col 16  p1               start on motif profile
col 17  p2               end on motif profile
col 18  q1               start on YOUR sequence (AA position)
col 19  q2               end on YOUR sequence (AA position)
col 20  env_from         envelope start
col 21  env_to           envelope end
col 22  acc              mean posterior probability of aligned residues
col 23  description      target sequence description
```

---

## 5. `Motif_psiblast_runner.sh`
**Purpose:** PSI-BLAST of motif consensus MSAs vs your RdRPs
**Argument style:** Positional (`$1 $2 ...`)

```bash
bash Motif_psiblast_runner.sh \
    $1  THREADS          \  # number of CPU threads
    $2  output_dir/      \  # output directory
    $3  input_fasta      \  # path to your_RdRPs.faa
    $4  Motif_type       \  # motif number: 1, 2, 3, or 4
    $5  motifDB_path     \  # path to msaFiles/ directory containing *Cons.*.faa files
    $6  out_prefix       \  # prefix for output .tsv filename
    $7  out_suffix       \  # suffix for output .tsv filename
    $8  rm_tmp              # "True" or "False" — remove temp files
```

**Output file:**
```
{out_prefix}_psiblast_{ini_name}_vs_mot{Motif_type}_{out_suffix}.tsv
```

**Example:**
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

**Output columns (PSI-BLAST outfmt 6 + profile name appended, 12 columns):**
```
col 1   sseqid           your sequence ID (subject)
col 2   pident           % identity
col 3   sstart           start on YOUR sequence = q1
col 4   send             end on YOUR sequence = q2
col 5   qstart           start on motif profile MSA = p1
col 6   qend             end on motif profile MSA = p2
col 7   slen             length of your sequence = qL
col 8   qlen             length of motif profile = pL
col 9   length           alignment length = ali_len
col 10  evalue           E-value
col 11  bitscore         bit score
col 12  profile_name     motif cluster name (appended from filename)
```

> Note: query/subject are **reversed** compared to other tools — the motif MSA is the query, your RdRP is the subject.

---

## 6. `Motif_MMseqs_runner.sh`
**Purpose:** MMseqs2 profile search of your RdRPs vs motif MMseqs2 profile DB
**Argument style:** Positional (`$1 $2 ...`)

```bash
bash Motif_MMseqs_runner.sh \
    $1  THREADS          \  # number of CPU threads
    $2  output_dir/      \  # output directory
    $3  input            \  # path to your_RdRPs.faa (or existing MMseqs2DB path)
    $4  Motif_type       \  # motif number: 1, 2, 3, or 4
    $5  motifDB_path     \  # path to MMseqs2 profile database (MMseqs2DB file, no extension)
    $6  out_prefix       \  # prefix for output .tsv filename
    $7  out_suffix       \  # suffix for output .tsv filename
    $8  rm_tmp              # "True" or "False" — remove temp files
```

**Output file:**
```
{out_prefix}_MMseqs2_vs_{ini_name}_vs_mot{Motif_type}_{out_suffix}.tsv
```

**Example:**
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

**Output columns (MMseqs2 convertalis custom format, 19 columns):**
```
col 1   query            your sequence ID
col 2   target           motif profile cluster name
col 3   evalue           E-value
col 4   gapopen          number of gap openings
col 5   pident           % identity
col 6   nident           number of identical residues
col 7   qstart           start on YOUR sequence = q1
col 8   qend             end on YOUR sequence = q2
col 9   qlen             length of your sequence = qL
col 10  tstart           start on motif profile = p1
col 11  tend             end on motif profile = p2
col 12  tlen             motif profile length = pL
col 13  alnlen           alignment length = ali_len
col 14  raw              raw score
col 15  bits             bit score
col 16  qframe           reading frame (0 for protein input)
col 17  mismatch         number of mismatches
col 18  qcov             query coverage
col 19  tcov             target (profile) coverage
```

---

## 7. `Motif_DiamondP_runner.sh`
**Purpose:** Diamond BLASTp — motif singletons as query vs your RdRPs as DB
**Argument style:** Positional (`$1 $2 ...`)

```bash
bash Motif_DiamondP_runner.sh \
    $1  THREADS          \  # number of CPU threads
    $2  output_dir/      \  # output directory
    $3  input            \  # path to your_RdRPs.faa (used as the SUBJECT database)
    $4  Motif_type       \  # motif number: 1, 2, 3, or 4
    $5  motifDB_path     \  # path to singletons.faa (motif singletons = QUERY)
    $6  out_prefix       \  # prefix for output .tsv filename
    $7  out_suffix       \  # suffix for output .tsv filename
    $8  rm_tmp              # "True" or "False" — remove temp files
```

**Output file:**
```
{out_prefix}_DiamondP_vs_{ini_name}_vs_mot{Motif_type}_{out_suffix}.tsv
```

**Example:**
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

**Output columns (Diamond BLASTp outfmt 6 custom, 14 columns):**
```
col 1   sseqid           your sequence ID (subject) = RdRp_id
col 2   qseqid           motif singleton ID, prefixed "mot.1." by awk = profile
col 3   evalue           E-value
col 4   gapopen          number of gap openings
col 5   pident           % identity
col 6   length           alignment length = ali_len
col 7   qstart           start on motif singleton (query) = p1
col 8   qend             end on motif singleton (query) = p2
col 9   sstart           start on YOUR sequence (subject) = q1
col 10  send             end on YOUR sequence (subject) = q2
col 11  qlen             length of motif singleton = pL
col 12  slen             length of your sequence = qL
col 13  mismatch         number of mismatches
col 14  bitscore         bit score
```

> Note: query/subject are **reversed** — motif singletons are query, your RdRPs are subject. `q1`/`q2` are in columns 9/10 (`sstart`/`send`), not 7/8.

---

## Quick reference: argument position summary

| Position | hhalign_runner (flag) | hhseach_runner (pos) | hmmsearch_runner (pos) | psiblast_runner (pos) | MMseqs_runner (pos) | DiamondP_runner (pos) |
|----------|-----------------------|----------------------|------------------------|----------------------|--------------------|-----------------------|
| Threads | `-t` | `$1` | `$1` | `$1` | `$1` | `$1` |
| Output dir | `-o` | `$2` | `$2` | `$2` | `$2` | `$2` |
| Input type | `-y` | `$3` | *(not used)* | *(not used)* | *(not used)* | *(not used)* |
| Input (your RdRPs) | `-i` | `$4` | `$3` | `$3` | `$3` | `$3` |
| Motif type (1–4) | `-M` | `$5` | `$4` | `$4` | `$4` | `$4` |
| MotifDB path | `-d` | `$6` | `$5` | `$5` | `$5` | `$5` |
| Output prefix | `-p` | `$7` | `$6` | `$6` | `$6` | `$6` |
| Output suffix | `-s` | `$8` | `$7` | `$7` | `$7` | `$7` |
| Remove tmps | `-r` | `$9` | `$8` | `$8` | `$8` | `$8` |
| E-value | `-e` | *(hardcoded 0.5)* | *(hardcoded 0.5)* | *(hardcoded 0.5)* | *(hardcoded 0.5)* | *(hardcoded 0.05)* |
| Extension | *(n/a)* | `$10` | *(n/a)* | *(n/a)* | *(n/a)* | *(n/a)* |

---

## MotifDB path per tool

Each tool expects a different type of motif database. These are all produced by `Profiler_motifs.sh`:

| Tool | `motifDB_path` argument points to |
|------|----------------------------------|
| `hhalign_runner.sh` | `mot.X/HHMfiles/db` (no file extension — HH-suite index prefix) |
| `hhseach_runner.sh` | `mot.X/HHMfiles/db` (same) |
| `hmmsearch_runner.sh` | `mot.X/HMMfiles/Vfin.mot.X_profiles.hmm` (HMMER3 `.hmm` file) |
| `psiblast_runner.sh` | `mot.X/msaFiles/` (directory containing `*Cons.*.faa` consensus MSA files) |
| `MMseqs_runner.sh` | `mot.X/MMseqs2_profiles/MMseqs2DB` (MMseqs2 DB prefix, no extension) |
| `DiamondP_runner.sh` | `mot.X/singletons.faa` (single merged FASTA of all singleton sequences) |

---

## Column name mapping to unified R table

After parsing with `GenericHitsParser()`, all tool outputs are standardised to these column names:

| Unified name | hhalign/hhseach | hmmsearch | psiblast | MMseqs2 | DiamondP |
|--------------|-----------------|-----------|----------|---------|---------|
| `q1` | col 7 | col 18 | col 3 (sstart) | col 7 (qstart) | col 9 (sstart) |
| `q2` | col 8 | col 19 | col 4 (send) | col 8 (qend) | col 10 (send) |
| `qL` | *(from faa)* | col 3 | col 7 (slen) | col 9 (qlen) | col 12 (slen) |
| `p1` | col 9 | col 16 | col 5 (qstart) | col 10 (tstart) | col 7 (qstart) |
| `p2` | col 10 | col 17 | col 6 (qend) | col 11 (tend) | col 8 (qend) |
| `pL` | col 4 | col 6 | col 8 (qlen) | col 12 (tlen) | col 11 (qlen) |
| `score` | col 13 | col 8 | col 11 (bitscore) | col 15 (bits) | col 14 (bitscore) |
| `evalue` | col 12 | col 7 | col 10 | col 3 | col 3 |
| `ali_len` | col 3 | *(q2−q1)* | col 9 (length) | col 13 (alnlen) | col 6 (length) |
| `subject_name` | col 1 | col 4 | col 12 (appended) | col 2 (target) | col 2 (qseqid) |
