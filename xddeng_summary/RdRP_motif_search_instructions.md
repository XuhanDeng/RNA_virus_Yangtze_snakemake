# How to Run `RdRP_motif_search.smk`

## What this pipeline does

Identifies RdRP motifs A, B, C, D in your candidate protein sequences using five parallel search tools (HHsearch, HMMER, PSI-BLAST, MMseqs2, Diamond). It replicates the RVMT motif identification approach without requiring any RVMT wrapper scripts.

---

## Prerequisites

### Software — must be in your `$PATH`

| Tool | Purpose | Install |
|------|---------|---------|
| `snakemake` | workflow engine | `conda install -c bioconda snakemake` |
| `hhconsensus`, `hhmake`, `hhsearch`, `splitfasta.pl`, `ffindex_build` | HH-suite | `conda install -c bioconda hhsuite` |
| `hmmbuild`, `hmmsearch` | HMMER3 | `conda install -c bioconda hmmer` |
| `mmseqs` | MMseqs2 | `conda install -c bioconda mmseqs2` |
| `makeblastdb`, `psiblast` | BLAST+ | `conda install -c bioconda blast` |
| `diamond` | Diamond BLASTp | `conda install -c bioconda diamond` |
| `parallel` | GNU Parallel | `conda install -c conda-forge parallel` |

Install all at once:
```bash
conda create -n rdRP_motif \
    snakemake hhsuite hmmer mmseqs2 blast diamond parallel \
    -c bioconda -c conda-forge -y
conda activate rdRP_motif
```

### Input files you must provide

| File | Description |
|------|-------------|
| `your_RdRPs.faa` | Your candidate RdRP protein sequences, multi-FASTA format |
| `Sequence_Library/mot.1–4/` | Pre-built MSA clusters from Zenodo (`RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library/`) |

---

## Setup

### 1. Edit `config.yaml`

Open [config.yaml](config.yaml) and set these four paths:

```yaml
input_faa: "/absolute/path/to/your_RdRPs.faa"
sequence_library: "/Users/dengxuhan/yangtze/software_test/RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library"
outdir: "/absolute/path/to/output_directory"
threads: 16        # set to your available CPU count
prefix: "MyRun"   # any short label you like
```

`sequence_library` already points to your Zenodo copy — only change it if the path is different.

### 2. Check your input FASTA

Sequences must be standard protein FASTA. IDs must not contain spaces:
```
>seq_001
MSLVKAQLLSAWKAHPTSTMT...
>seq_002
ELNTLVNRGYGAVNWQTEREH...
```

---

## Run

```bash
# Activate environment
conda activate rdRP_motif

# Dry run first — shows all jobs without executing
snakemake -s RdRP_motif_search.smk --configfile config.yaml -j 16 -n

# Full run
snakemake -s RdRP_motif_search.smk --configfile config.yaml -j 16
```

`-j 16` = maximum 16 jobs in parallel. Set to match your CPU count.

---

## Pipeline steps and output files

```
Sequence_Library/mot.1–4/*.afa
        │
        ▼ rule add_consensus  (×4, one per motif)
        │   hhconsensus -M 50 -cov 50
        │   → outdir/profiles/mot.{1-4}/msaFiles/*.Cons.msa.afa
        │
        ├──▶ rule build_hmm          (×4)
        │     hmmbuild → cat *.hmm
        │     → outdir/profiles/mot.{1-4}/HMMfiles/profiles.hmm
        │          └───────────────▶ rule search_hmmsearch  (×4)
        │                             → outdir/search/hmmsearch/mot.{1-4}.tsv
        │
        ├──▶ rule build_hhm          (×4)
        │     hhmake → ffindex_build
        │     → outdir/profiles/mot.{1-4}/HHMfiles/db
        │          ├───────────────▶ rule search_hhalign    (×4)
        │          │                  → outdir/search/hhalign/mot.{1-4}.tsv
        │          └──▶ rule build_mmseqs_profile  (×4)
        │                mmseqs convertprofiledb
        │                → outdir/profiles/mot.{1-4}/MMseqs2_profiles/MM
        │                     └────▶ rule search_mmseqs     (×4)
        │                             → outdir/search/mmseqs/mot.{1-4}.tsv
        │
        ├──▶ rule collect_singletons (×4)
        │     → outdir/profiles/mot.{1-4}/singletons.faa
        │          └───────────────▶ rule search_diamond    (×4)
        │                             → outdir/search/diamond/mot.{1-4}.tsv
        │
        └──────────────────────────▶ rule search_psiblast   (×4)
                                      → outdir/search/psiblast/mot.{1-4}.tsv
```

### Final output structure

```
outdir/
├── profiles/
│   ├── mot.1/
│   │   ├── msaFiles/          # *.Cons.msa.afa — consensus MSAs
│   │   ├── HMMfiles/
│   │   │   └── profiles.hmm   # HMMER3 profile DB for motif A
│   │   ├── HHMfiles/
│   │   │   ├── db             # HH-suite profile DB (ffindex data)
│   │   │   └── db.index       # HH-suite profile DB (ffindex index)
│   │   ├── MMseqs2_profiles/
│   │   │   └── MM             # MMseqs2 profile DB for motif A
│   │   └── singletons.faa     # single-sequence clusters (for Diamond)
│   ├── mot.2/  ...
│   ├── mot.3/  ...
│   └── mot.4/  ...
├── search/
│   ├── hhalign/
│   │   ├── mot.1.tsv          # 13-column HHsearch blasttab output
│   │   ├── mot.2.tsv
│   │   ├── mot.3.tsv
│   │   └── mot.4.tsv
│   ├── hmmsearch/
│   │   └── mot.{1-4}.tsv      # 23-column HMMER domtblout (comments stripped)
│   ├── psiblast/
│   │   └── mot.{1-4}.tsv      # 12-column PSI-BLAST outfmt6 + cluster name
│   ├── mmseqs/
│   │   └── mot.{1-4}.tsv      # 19-column MMseqs2 convertalis output
│   └── diamond/
│       └── mot.{1-4}.tsv      # 14-column Diamond outfmt6, col2 prefixed "mot.N."
└── logs/
    └── *.log                   # per-rule log files for debugging
```

---

## Output column reference

> **Key columns shared across all tools after R parsing:**
> `q1`/`q2` = start/end of the motif on **your** sequence (AA position)
> `p1`/`p2` = start/end on the motif profile
> `score`, `evalue` = hit quality

| Tool | File | `q1` col | `q2` col | `score` col | `evalue` col |
|------|------|----------|----------|-------------|-------------|
| HHsearch | `hhalign/mot.N.tsv` | 7 | 8 | 13 | 12 |
| HMMER | `hmmsearch/mot.N.tsv` | 18 | 19 | 8 | 7 |
| PSI-BLAST | `psiblast/mot.N.tsv` | 3 (sstart) | 4 (send) | 11 | 10 |
| MMseqs2 | `mmseqs/mot.N.tsv` | 7 (qstart) | 8 (qend) | 15 (bits) | 3 |
| Diamond | `diamond/mot.N.tsv` | 9 (sstart) | 10 (send) | 14 | 3 |

PSI-BLAST and Diamond have **reversed query/subject** — the motif sequence is the query and your RdRP is the subject, so `q1`/`q2` are in the subject coordinate columns (`sstart`/`send`).

---

## Next step: parse in R

Load all 20 output files into R using `GenericHitsParser()` from `Utils/basicf.r`:

```r
source("Utils/basicf.r")

tools  <- c("hhalign", "hmmsearch", "psiblast", "mmseqs", "diamond")
motifs <- 1:4

all_hits <- rbindlist(lapply(motifs, function(m) {
    rbindlist(lapply(tools, function(tool) {
        f <- sprintf("outdir/search/%s/mot.%d.tsv", tool, m)
        h <- GenericHitsParser(
            inpt          = f,
            search_tool   = tool,
            input_was     = "new_name",
            Query2Profile = TRUE,
            breakhdrs     = FALSE
        )
        h$motif_type <- m
        h
    }), fill = TRUE)
}), fill = TRUE)

# Apply Neri's final filter thresholds
all_hits <- all_hits[evalue < 1e-5]
all_hits <- all_hits[score  > 10  ]
all_hits <- all_hits[abs(q2 - q1) > 10]   # ali_len on query > 10 AA
```

See [Neri_filter_standards.md](Neri_filter_standards.md) for the complete filtering rationale and [Motif_identification_your_RdRPs.md](Motif_identification_your_RdRPs.md) for the full R parsing workflow.

---

## Troubleshooting

**`splitfasta.pl` not found**
It ships with HH-suite. Check: `which splitfasta.pl`. If missing, reinstall hhsuite or add `$CONDA_PREFIX/bin` to your PATH.

**`ffindex_build` not found**
Also part of HH-suite. Check: `conda install -c bioconda hhsuite`.

**MMseqs2 out of memory**
Reduce `--split-memory-limit` in the `search_mmseqs` rule (line ~513) from `100G` to match your available RAM.

**PSI-BLAST slow**
The psiblast rule loops over ~300 cluster files per motif sequentially. This is inherently slow. Run it overnight or reduce the Sequence_Library to only the most relevant clusters.

**Rerun a single rule after fixing an error**
```bash
# Force rerun of one specific rule for motif 2
snakemake -s RdRP_motif_search.smk --configfile config.yaml \
    -j 16 --force-rerun search_psiblast \
    outdir/search/psiblast/mot.2.tsv
```

**Visualise the DAG**
```bash
snakemake -s RdRP_motif_search.smk --configfile config.yaml \
    --dag | dot -Tpdf > dag.pdf
```
