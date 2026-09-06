# Phylogenetic Tree Construction

Goal: insert newly identified RdRp sequences (from 48 RNA samples + ICTV references) into
the Neri et al. phylogenetic framework. Only palm-confirmed categories are used because
palm_annot provides A/B/C motif positions in FASTA headers, which are required for
de-permutation.

---

## Input categories (palm-confirmed only)

| Category | Source |
|---|---|
| `High-confident` | RC + LP + palm confirmed |
| `RDRPCatch_palm` | RC + palm confirmed |
| `lucaprot_palmscan` | LP + palm confirmed |

---

## Step 1 — Collect input sequences

Copy the per-category protein FASTAs into the input folder without modification.

**Contig sequences** → `result/07_phylogenetic_tree/0_input/Contig/`
```
result/03_RDRP_identification/8_RdRp_protein/High-confident.faa
result/03_RDRP_identification/8_RdRp_protein/RDRPCatch_palm.faa
result/03_RDRP_identification/8_RdRp_protein/lucaprot_palmscan.faa
```

**ICTV sequences** → `result/07_phylogenetic_tree/0_input/ICTV/`
```
/scratch/…/result/03_RDRP_identification/ICTV/8_RdRp_protein/High-confident.faa
/scratch/…/result/03_RDRP_identification/ICTV/8_RdRp_protein/RDRPCatch_palm.faa
/scratch/…/result/03_RDRP_identification/ICTV/8_RdRp_protein/lucaprot_palmscan.faa
```

Script: `scripts/07_phylogenetic_tree/01_collect_input.py` (or simple cp/cat in snakemake rule)

---

## Step 2 — De-permute RdRp sequences

Some RdRps have a permuted motif order (C-A-B instead of canonical A-B-C). These must be
re-arranged before alignment.

### Header format (from palm_annot)

```
>{sequence_id} A:{pos}:{motif} B:{pos}:{motif} C:{pos}:{motif}
```

Example:
```
>GouBa_N_1_0000000057_frame=1_RdRp_681-977 A:199:IASDVSGWEKNF B:267:SGCLMTTSSNGVAR C:96:YEGDFSEF
```

Positions are 1-based AA indices.

### Detection

Parse `A_pos`, `B_pos`, `C_pos` from header.

- `C_pos < A_pos` → **permuted** (C-A-B order), needs de-permutation
- `A_pos < B_pos < C_pos` → **canonical** (A-B-C order), write as-is

### De-permutation logic

Palm_annot reports **1-based** positions. Convert to 0-based before slicing:

```python
A_pos = int(A_pos_str) - 1
```

Cut at `A_pos`. Everything before A (the C-containing N-terminal block) moves to the end:

```
original:   [0 → A_pos)        [A_pos → end)
            C block             A-B block

de-permuted: [A_pos → end) + [0 → A_pos)
             A-B block         C block
```

```python
seq = "".join(line.strip() for line in raw_lines)   # clean: no whitespace/tabs
seq_ab = seq[A_pos:]
seq_c  = seq[:A_pos]
depermuted = seq_ab + seq_c
```

Verify motifs are present after de-permutation (sanity check):
```python
assert A_motif in depermuted
assert B_motif in depermuted
assert C_motif in depermuted
```

### Output format

- Header: simplified to `>{sequence_id}` only (drop A/B/C annotation)
- Sequence: reformatted to 60 AA per line, no tabs

```
>GouBa_N_1_0000000057_frame=1_RdRp_681-977
IKFPAGAIEAAQQGVTEMYKAAGFVWRCGFVGTEHGEALARRFDEVYPTIV
EDVCKHSKPGYPYRLAYGTNEKLLQERPDWVRQLVWQRVKNIAEKADRKLW
...
```

### Output paths

```
result/07_phylogenetic_tree/1_depermuted/Contig/   — per-category depermuted FASTAs
result/07_phylogenetic_tree/1_depermuted/ICTV/     — per-category depermuted FASTAs
```

### Merge

After de-permutation, concatenate all Contig + ICTV FASTAs into one file:

```
result/07_phylogenetic_tree/1_depermuted/Contig_ICTV_rdrp.faa
```

### Summary table

After merging, write a TSV summary table:

```
result/07_phylogenetic_tree/1_depermuted/Contig_ICTV_rdrp_summary.tsv
```

One row per sequence. Columns:

| Column | Content |
|---|---|
| `seq_id` | sequence ID (simplified header, no A/B/C annotation) |
| `source` | `Contig` or `ICTV` |
| `category` | `High-confident`, `RDRPCatch_palm`, or `lucaprot_palmscan` |
| `A_pos` | A motif start position (1-based, from palm_annot) |
| `A_motif` | A motif amino acid sequence |
| `B_pos` | B motif start position (1-based) |
| `B_motif` | B motif amino acid sequence |
| `C_pos` | C motif start position (1-based) |
| `C_motif` | C motif amino acid sequence |
| `depermuted` | `yes` if C_pos < A_pos (was permuted and fixed), `no` if canonical |

Example rows:
```
seq_id                                          source  category         A_pos  A_motif       B_pos  B_motif         C_pos  C_motif   depermuted
GouBa_N_1_0000000057_frame=1_RdRp_681-977      Contig  RDRPCatch_palm   199    IASDVSGWEKNF  267    SGCLMTTSSNGVAR  96     YEGDFSEF  yes
NC_004102_frame=1_RdRp_1-300                   ICTV    High-confident   45     FADDMSR       112    TGLLSSSGVAR     189    YEGDFSE   no
```

Script: `scripts/07_phylogenetic_tree/02_depermute.py`

---

## Step 3 — Add Neri et al. reference sequences

Merge Yangtze + ICTV sequences with the RVMT reference RdRp dataset used in
Neri et al. (2022) for placement into the existing phylogenetic framework.

### Reference dataset

RVMT database RdRp sequences (already downloaded):
```
database/RVMT/RVMT_Zenodo_V4/
```

Extract the RdRp palm-core sequences from RVMT. The RVMT FASTA contains full-length
protein predictions — use the palm-trimmed version if available, otherwise use full-length.

```
database/RVMT/RVMT_Zenodo_V4/RVMT_RdRp.faa   ← source (check actual filename in dir)
```

Output:
```
result/07_phylogenetic_tree/2_with_ref/Contig_ICTV_RVMT_rdrp.faa
```

This file is the union of:
- `result/07_phylogenetic_tree/1_depermuted/Contig_ICTV_rdrp.faa`
- RVMT RdRp sequences (headers prefixed with `RVMT|` for easy identification)

Script: `scripts/07_phylogenetic_tree/03_add_reference.py`

---

## Step 4 — Multiple sequence alignment (MAFFT)

Align the merged sequence set with MAFFT.

```
result/07_phylogenetic_tree/2_with_ref/Contig_ICTV_RVMT_rdrp.faa
  → result/07_phylogenetic_tree/3_aligned/Contig_ICTV_RVMT_rdrp_aligned.faa
```

Command:
```bash
mafft --auto --thread {threads} \
    result/07_phylogenetic_tree/2_with_ref/Contig_ICTV_RVMT_rdrp.faa \
    > result/07_phylogenetic_tree/3_aligned/Contig_ICTV_RVMT_rdrp_aligned.faa
```

`--auto` selects L-INS-i for < ~200 sequences, FFT-NS-2 for larger datasets.

---

## Step 5 — Alignment trimming (TrimAl)

Remove poorly aligned columns.

```
result/07_phylogenetic_tree/3_aligned/Contig_ICTV_RVMT_rdrp_aligned.faa
  → result/07_phylogenetic_tree/4_trimmed/Contig_ICTV_RVMT_rdrp_trimmed.faa
```

Command:
```bash
trimal -in  result/07_phylogenetic_tree/3_aligned/Contig_ICTV_RVMT_rdrp_aligned.faa \
       -out result/07_phylogenetic_tree/4_trimmed/Contig_ICTV_RVMT_rdrp_trimmed.faa \
       -automated1
```

`-automated1`: heuristic that selects `-gappyout` or `-strict` depending on alignment
length, recommended for phylogenetics.

---

## Step 6 — Phylogenetic tree (IQ-TREE2)

Build maximum likelihood tree with automatic model selection.

```
result/07_phylogenetic_tree/4_trimmed/Contig_ICTV_RVMT_rdrp_trimmed.faa
  → result/07_phylogenetic_tree/5_tree/Contig_ICTV_RVMT_rdrp.treefile
  → result/07_phylogenetic_tree/5_tree/Contig_ICTV_RVMT_rdrp.iqtree   (log + model info)
```

Command:
```bash
iqtree2 -s result/07_phylogenetic_tree/4_trimmed/Contig_ICTV_RVMT_rdrp_trimmed.faa \
        -m TEST \
        -B 1000 \
        -T {threads} \
        --prefix result/07_phylogenetic_tree/5_tree/Contig_ICTV_RVMT_rdrp
```

`-m TEST`: ModelTest-NG model selection  
`-B 1000`: 1000 ultrafast bootstrap replicates

---

## Output structure

```
result/07_phylogenetic_tree/
  0_input/
    Contig/
      High-confident.faa
      RDRPCatch_palm.faa
      lucaprot_palmscan.faa
    ICTV/
      High-confident.faa
      RDRPCatch_palm.faa
      lucaprot_palmscan.faa

  1_depermuted/
    Contig/
      High-confident.faa
      RDRPCatch_palm.faa
      lucaprot_palmscan.faa
    ICTV/
      High-confident.faa
      RDRPCatch_palm.faa
      lucaprot_palmscan.faa
    Contig_ICTV_rdrp.faa              ← merged, de-permuted
    Contig_ICTV_rdrp_summary.tsv      ← per-sequence motif + permutation table

  2_with_ref/
    Contig_ICTV_RVMT_rdrp.faa         ← merged with RVMT references

  3_aligned/
    Contig_ICTV_RVMT_rdrp_aligned.faa

  4_trimmed/
    Contig_ICTV_RVMT_rdrp_trimmed.faa

  5_tree/
    Contig_ICTV_RVMT_rdrp.treefile
    Contig_ICTV_RVMT_rdrp.iqtree
    Contig_ICTV_RVMT_rdrp.log
```

---

## Scripts

| Step | Script |
|---|---|
| 1 | `scripts/07_phylogenetic_tree/01_collect_input.py` |
| 2 | `scripts/07_phylogenetic_tree/02_depermute.py` ← **done** |
| 3 | `scripts/07_phylogenetic_tree/03_add_reference.py` |
| 4–6 | Snakemake rules calling mafft / trimal / iqtree2 directly |

---

## Example

### Permuted input (C-A-B order, C_pos=96 < A_pos=199 < B_pos=267):
```
>GouBa_N_1_0000000057_frame=1_RdRp_681-977 A:199:IASDVSGWEKNF B:267:SGCLMTTSSNGVAR C:96:YEGDFSEF
IKFPAGAIEAAQQGVTEMYKAAGFVWRCGFVGTEHGEALARRFDEVYPTIVEDVCKHSKPGYPYRLAYGTNEKLLQERPD
WVRQLVWQRVKNIAEYEGDFSEFKADRKLWILRDMRDPVRLFGKNQGQPIRKPMCRIISHVSLIDQMVMRFFFGAYAAAE
GEFYPYLPTKKGIGFSSEHATKIGNCVYATSQELDRDPIASDVSGWEKNFSQDCAEIFARHMLATCEDRCSLLEKAASWW
KESLTSTPYVTDGGQLIDYADTRVQRSGCLMTTSSNGVARVACAVACDHIANAMGDD
```

### De-permuted output (canonical A-B-C order, header simplified):

A_pos = 199 - 1 = 198 (0-based). Cut: `seq[198:] + seq[:198]`

```
>GouBa_N_1_0000000057_frame=1_RdRp_681-977
IASDVSGWEKNFSQDCAEIFARHMLATCEDRCSLLEKAASWWKESLTSTPYVTDGGQLIDYADTRVQRSGCLMTTSSNGV
ARVACAVACDHIANAMGDDIKFPAGAIEAAQQGVTEMYKAAGFVWRCGFVGTEHGEALARRFDEVYPTIVEDVCKHSKPG
YPYRLAYGTNEKLLQERPDWVRQLVWQRVKNIAEYEGDFSEFKADRKLWILRDMRDPVRLFGKNQGQPIRKPMCRIISHV
SLIDQMVMRFFFGAYAAAEGEFYPYLPTKKGIGFSSEHATKIGNCVYATSQELDRDP
```

---

## Notes

- All three input categories (`High-confident`, `RDRPCatch_palm`, `lucaprot_palmscan`) have
  A/B/C positions directly in the FASTA header from palm_annot. Header format is identical:
  `>{sequence_id} A:{pos}:{motif} B:{pos}:{motif} C:{pos}:{motif}`
- `lucaprot_palmscan` headers use ORF-style sequence IDs (`ORF{n}_{contig}:{start}:{end}`)
  but still carry full A/B/C annotation — no TSV lookup needed.
- RVMT sequences are already in canonical A-B-C order; no de-permutation needed for them.
- Sequence IDs must be unique across Contig + ICTV + RVMT. Contig IDs contain sample name
  (e.g. `GouBa_N_1_…`), ICTV IDs contain accession (e.g. `NC_004102_…`), RVMT IDs are
  prefixed with `RVMT|` to avoid collision.
