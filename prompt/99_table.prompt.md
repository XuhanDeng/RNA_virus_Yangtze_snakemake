
## esvirtu_info_table
2_esvirtu_table

### Input files

**input_table1**: `database/esviritu/v3.2.4/virus_pathogen_database.all_metadata.tsv`
columns: Accession, description, Name, Segment, kingdom, phylum, tclass, order, family, genus, species, subspecies, Length, TaxID, Assembly, Asm_length

**input_table2 (RPKMF)**: `result/99_tables/01_cp_table/04_modified_esvirtue_ribodetector/esvirtu_rpkmf.subspecies.tsv`
columns: subspecies (index), then 48 sample RPKMF columns

**input_table2b (read count)**: `result/99_tables/01_cp_table/04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv`
columns: subspecies (index), then 48 sample read_count columns, last column: Asm_length

**input_table3 (virushostdb)**: `database/virushostdb/virushostdb.tsv`
columns: virus tax id, virus name, virus lineage, refseq id, KEGG GENOME, KEGG DISEASE, DISEASE, host tax id, host name, host lineage, pmid, evidence, sample type, source organism

**input_table4 (host_metadata)**: `database/esviritu/v3.2.4/virus_pathogen_database.host_metadata.v3.2.4.tsv`
columns: Accession, description, Name, Segment, kingdom, phylum, tclass, order, family, genus, species, subspecies, Length, TaxID, Assembly, Asm_length, Host

**output**: `result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv`

---

### Operations

**Step 0 — filter input_table2 and input_table2b**
Keep only rows where subspecies starts with `t__` (ESviritu-identified). Remove `unknown::acc:` assembly contigs.
Extract `Asm_length` from input_table2b; it will be added as a column in the output.

**Step 1 — join taxonomy from input_table1**
Use subspecies column as join key. Extract columns: Name, Segment, kingdom, phylum, tclass, order, family, genus, species, TaxID.
- One subspecies may match multiple rows in input_table1; values should be consistent — warn if not, take first.
- Accession: one subspecies may have multiple accessions → join all into one cell separated by `;`.

**Step 2 — add genome_type from phylum**
Map phylum → genome_type using ICTV classification:
- ssDNA: Cressdnaviricota, Cossaviricota, Preplasmiviricota
- ssRNA(+): Pisuviricota, Kitrinoviricota, Lenarviricota
- ssRNA(-): Negarnaviricota
- dsRNA: Duplornaviricota
- dsDNA: Nucleocytoviricota, Peploviricota, Uroviricota, Dividoviricota, Helvetiaviricota
- dsDNA-RT: Artverviricota
- unknown: unclassified_Viruses or unmapped
Print summary count per genome_type.

**Step 3 — join host from virushostdb (input_table3)**
Join on TaxID (virus tax id). Add columns: vhd_virus_name, vhd_host_tax_id, vhd_host_name, vhd_host_lineage, vhd_pmid, vhd_evidence.
- Exclude rows where host tax id == 1 or host name == "root" before aggregating.
- Multiple host entries per TaxID → join with `;`.

**Step 4 — join host from host_metadata (input_table4)**
Join on first Accession (split `;`, take first). Add column: host_meta_Host.
- Exclude rows where Host is "root", "NA", or empty.

**Step 5 — merged_host**
Add column `merged_host`: use vhd_host_name if non-empty, else fall back to host_meta_Host.

**Step 6 — column order**
Final column order:
subspecies | Name | Segment | kingdom | phylum | tclass | order | family | genus | species | TaxID | genome_type | Asm_length | Accession | vhd_virus_name | vhd_host_tax_id | vhd_host_name | vhd_host_lineage | vhd_pmid | vhd_evidence | host_meta_Host | merged_host | [48 sample RPKMF columns] | [48 sample read_count columns]

### Script
`scripts/99_tables/esvirtu_info_table.py`

Run:
```bash
python scripts/99_tables/esvirtu_info_table.py
```




## 3_RNA_virus_table

### Input files

**input_table1 (contig abundance, RPKMF)**:
`result/99_tables/01_cp_table/04_modified_esvirtue_ribodetector/esvirtu_rpkmf.subspecies.tsv`
- Use only rows where first column starts with `unknown::acc:` (assembly contigs, not ESvirtu references)
- Strip `unknown::acc:` prefix to get bare contig ID (e.g. `AnQing_D_33_0000000012`)
- columns: contig_id (derived), then 48 sample RPKMF columns

**input_table1b (contig abundance, read count)**:
`result/99_tables/01_cp_table/04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv`
- Same row filter: keep only `unknown::acc:` rows, strip prefix to get contig_id
- columns: contig_id (derived), then 48 sample read_count columns, last column: Asm_length (contig length)

**input_table2 (RdRp identification)**:
`result/99_tables/01_cp_table/03_RDRP_identification/tier_summary.tsv`
- key column for join: `contig`
- columns to add: `RdRp_id`, `tier`, `motif_order`, `n_motifs`
- If one contig has multiple rows, keep the highest tier: tier1a > tier1b > tier2 > tier3

**input_table3 (bait recovery)**:
`result/99_tables/01_cp_table/06_contig_bait/bait_target_map.tsv`
- key column for join: `recovered_contig`
- columns: recovered_contig, bait_source, bait_category, best_bait_contig, pident, tcov, evalue, n_bait_hits
- columns to add: `bait_category`, `best_bait_contig`, `bait_source`, `pident`, `tcov`, `evalue`, `n_bait_hits`
- If one contig has multiple bait rows, keep highest confidence: tier1a > tier1b > tier2 > RVMT

**output**: `result/99_tables/3_RNA_virus_table/rna_virus_contig_table.tsv`

---

### Operations

**Step 1 — load contig abundance (input_table1 + input_table1b)**
Filter rows starting with `unknown::acc:` from both tables. Strip prefix to get `contig_id`.
From input_table1b, extract `Asm_length` (last column = contig length).

**Step 2 — join RdRp info (input_table2)**
Join on `contig_id` == `contig`. Add columns: `RdRp_id`, `tier`, `motif_order`, `n_motifs`.
If multiple RdRp rows match one contig, keep the row with highest tier (tier1a > tier1b > tier2 > tier3). Contigs with no match: leave these columns empty.

**Step 3 — join bait info (input_table3)**
Join on `contig_id` == `recovered_contig`. Add columns: `bait_category`, `best_bait_contig`, `bait_source`, `pident`, `tcov`, `evalue`, `n_bait_hits`.
If multiple bait rows match one contig, keep highest confidence (tier1a > tier1b > tier2 > RVMT). Contigs with no match: leave empty.

**Step 4 — add contig_tag column**
Assign a single tag per contig based on combined RdRp + bait evidence, ranked:
`tier1a` > `tier1b` > `tier2` > `bait_tier1a` > `bait_tier1b` > `bait_tier2` > `bait_RVMT` > `tier3` > `unknown`

Rules:
- If RdRp tier is tier1a/tier1b/tier2 → use that tier directly (top priority)
- Else if bait_category is tier1a/tier1b/tier2/RVMT → use `bait_<category>` (bait ranks above tier3)
- Else if RdRp tier is tier3 → use `tier3`
- Else → `unknown`

**Step 5 — column order**
contig_id | Asm_length | RdRp_id | tier | motif_order | n_motifs | bait_category | best_bait_contig | bait_source | pident | tcov | evalue | n_bait_hits | contig_tag | [48 sample RPKMF columns] | [48 sample read_count columns last]

### Script
`scripts/99_tables/rna_virus_contig_table.py`

---


## 4_SSDNA_virus_table
will be supplenmented


## 5_Recalculated_RPKM_result

All outputs in this section are saved under `result/99_tables/5_Recalculated_RPKM_result/`.

Goal: Recalculate RPKMF using only RNA virus sequences as the reference set (instead of all detected sequences). This gives a more accurate relative abundance because non-RNA-virus sequences are excluded from the denominator.

---

### 1_RNA_virus_Recalculated_RPKM

**Objective**: Recalculate RPKMF restricted to RNA virus sequences only. The two reference sets are combined first into a single table, then RPKMF is calculated using the combined total as the denominator (so all values are comparable on a shared scale).

The reference set consists of:
1. **ESvirtu RNA virus rows** — rows in `esvirtu_info_table.tsv` where `genome_type` is `ssRNA(+)`, `ssRNA(-)`, or `dsRNA`.
2. **Novel RNA virus contigs** — rows in `rna_virus_contig_table.tsv` where `contig_tag` exactly matches one of (use exact string match; do not use a catch-all — new tags may be added later):
   - `tier1a`, `tier1b`, `tier2`, `tier3`
   - `bait_tier1a`, `bait_tier1b`, `bait_tier2`, `bait_RVMT`

---

**Input table1 (ESvirtu table)**
`result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv`
- Filter: keep rows where `genome_type` in `{ssRNA(+), ssRNA(-), dsRNA}`
- All columns are kept. Row ID column: `subspecies`.
- Metadata columns: `subspecies`, `Name`, `Segment`, `kingdom`, `phylum`, `tclass`, `order`, `family`, `genus`, `species`, `TaxID`, `genome_type`, `Asm_length`, `Accession`, `vhd_virus_name`, `vhd_host_tax_id`, `vhd_host_name`, `vhd_host_lineage`, `vhd_pmid`, `vhd_evidence`, `host_meta_Host`, `merged_host`
- Then 48 sample `*_read_count` columns (used for RPKMF calculation)
- Columns not present in table2 will be filled with `NA` in the combined table.

**Input table2 (RNA virus contig table)**
`result/99_tables/3_RNA_virus_table/rna_virus_contig_table.tsv`
- Filter: keep rows where `contig_tag` in the exact set listed above
- All columns are kept. Row ID column: `contig_id`.
- Metadata columns: `contig_id`, `Asm_length`, `RdRp_id`, `tier`, `motif_order`, `n_motifs`, `bait_category`, `best_bait_contig`, `bait_source`, `pident`, `tcov`, `evalue`, `n_bait_hits`, `contig_tag`
- Then 48 sample `*_read_count` columns (used for RPKMF calculation)
- Columns not present in table1 will be filled with `NA` in the combined table.

---

**Operations**

**Step 1 — combine**
Stack input table1 (filtered) and input table2 (filtered) into a single table.
- Rename `subspecies` / `contig_id` to a common column `seq_id`.
- Add a `source` column: `esvirtu` for rows from table1, `contig` for rows from table2.
- All metadata columns from both tables are retained; columns absent in one source are filled with `NA`.
- Both tables share the same 48 sample `*_read_count` column names — align on column name.

**Step 2 — recalculate RPKMF**
For each sample column `S`:
```
RPKMF(seq_i, S) = (read_count(seq_i, S) / Asm_length(seq_i) / 1000) /
                  (sum over all seq_j of read_count(seq_j, S) / 1e6)
```
"Million mapped reads" is the total mapped to the combined RNA virus reference set only (not all viruses).

**Step 3 — output**
Write two output files:

- `result/99_tables/5_Recalculated_RPKM_result/1_RNA_virus_Recalculated_RPKM/rna_virus_recalc_rpkmf.tsv`
  Columns: `seq_id | source | [all metadata cols from both tables] | Asm_length | [48 sample recalculated RPKMF columns]`
  Sample column names: replace `_read_count` → `_rpkmf`.

- `result/99_tables/5_Recalculated_RPKM_result/1_RNA_virus_Recalculated_RPKM/rna_virus_recalc_read_count.tsv`
  Same structure but with raw read counts instead of RPKMF.
  Columns: `seq_id | source | [all metadata cols from both tables] | Asm_length | [48 sample read_count columns]`

Print per-sample total mapped reads (sum of read counts across all RNA virus rows) as a QC check.

**Script**
`scripts/99_tables/rna_virus_recalc_rpkmf.py`


---

## 6_RNA_virus_Recalculated_correlation

Run the same Spearman all-vs-all correlation pipeline as workflow 05, but using the **recalculated RPKMF table** (`rna_virus_recalc_rpkmf.tsv`) as input instead of the original ESvirtu RPKMF. This restricts the correlation to RNA virus sequences only (both ESvirtu references and novel contigs combined).

All outputs saved under `result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation/`.
All logs under `log/99_tables/6_RNA_virus_Recalculated_correlation/`.

Reuse existing scripts from `scripts/05_esviritu/` — no new scripts needed.

---

### Input

**Recalculated RPKMF table**:
`result/99_tables/5_Recalculated_RPKM_result/1_RNA_virus_Recalculated_RPKM/rna_virus_recalc_rpkmf.tsv`
- Row ID column: `seq_id` (mix of ESvirtu subspecies names `t__...` and contig IDs)
- 48 sample `*_rpkmf` columns

**RdRp tier summary** (for annotation step):
`result/03_RDRP_identification/4_motif_search/motif_results/tier_summary.tsv`

---

### Steps

**Step 1 — Spearman all-vs-all** (reuse `scripts/05_esviritu/spearman_all_vs_all.py`)
- Input: recalculated RPKMF table
- Note: the row ID column is `seq_id`, not ending in `_RPKMF` — the script uses the first non-`_RPKMF` column as the label column, so this works as-is
- Output: `result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation/2_all_vs_all/rna_virus.all_vs_all.spearman.filtered.tsv`
- Use same config parameters as workflow 05: `min_samples`, `thresholds`, `p_threshold`, `chunk_size_all_vs_all`, `threads`

**Step 2 — filter into known/unknown pairs** (reuse `scripts/05_esviritu/filter_spearman_pairs.py`)
- Input: filtered spearman TSV from step 1
- Known: ESvirtu rows (`t__` prefix); Unknown: contig rows (bare contig IDs, no `t__` prefix)
- Note: the existing script classifies unknown by `unknown::` or `acc:` prefix — contig IDs in this table are bare (no prefix). The script may need a small adjustment or a wrapper to correctly classify `t__` rows as known and everything else as unknown.
- Output directory: `result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation/3_spearman_analysis/`
- Thresholds: 0.6, 0.7, 0.8, 0.9
- Pair types: `known_known_pair`, `known_unknown_pair`

**Step 3 — annotate known-unknown pairs with RdRp tier** (reuse `scripts/05_esviritu/annotate_rdrp_category.py`)
- Input: `known_unknown_pair/r{threshold}.tsv` from step 2 + `tier_summary.tsv`
- Output: `result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation/4_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv`
- Thresholds: 0.6, 0.7, 0.8, 0.9

---

### Snakemake workflow

Add rules to `workflow/99_tables.smk` (do not create a separate smk file).
Thresholds: `["0.6", "0.7", "0.8", "0.9"]`
Pair types: `["known_known_pair", "known_unknown_pair"]`
