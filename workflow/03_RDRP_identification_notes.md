# 03_RDRP_identification — redesign notes

## Why this changed

The old pipeline classified every contig into one of 8 tiers (`High-confident`,
`RDRPCatch_palm`, `RDRPCatch_lucaprot`, `RDRPCatch_only`, `lucaprot_palmscan`,
`lucaprot_only`, plus unions) based on whether RdRpCATCH and/or LucaProt hit it,
and whether palmscan separately confirmed a palmprint motif on each source.
It also used an overlap/strand check to decide whether an LP ORF and an RC
region were "the same" RdRp region on a contig.

New requirement: drop the tier system entirely. Keep only regions where
palmscan confirms an **exact ABC or CAB** palmprint motif. Still use overlap
to deduplicate (RdRpCATCH wins ties), but instead of merging overlapping
regions into a "same_region" super-category, just keep the single **longest**
surviving region per contig. Later requirement: don't lose the fact that
LucaProt *also* identified a contig even when its region lost the overlap
tie-break to RdRpCATCH — needed to evaluate LucaProt's identification
ability independently of which source ends up "winning" per contig.

## New logic (`merge_rdrp_results.py`)

1. **Load RdRpCATCH regions** — one row per contig (RC gives at most one
   region). Region coordinates (`aa_start`, `aa_end`, `frame`) come straight
   from RdRpCATCH's own output columns (`RdRp_from(AA)`, `RdRp_to(AA)`,
   `Translated_sequence_name (frame)`) — **not** from palmscan.
2. **Load LucaProt regions** — possibly multiple ORFs per contig. Coordinates
   (`orf_start`, `orf_end`) come from parsing LucaProt's `protein_id` field
   (nt coordinates), again not from palmscan.
3. **Motif filter** — for both sources, join in palmscan's own output
   (`pssm_ABC` column) as a pass/fail check only. A region survives only if
   `pssm_ABC` is exactly `ABC` or `CAB` (case-insensitive). Palmscan never
   supplies coordinates — it only confirms/rejects a candidate region that
   RdRpCATCH or LucaProt already defined.
4. **Dedup by overlap, RC wins, but nothing is discarded from the record** —
   for each contig, any motif-confirmed LP ORF that overlaps (same strand,
   ≥50% of the RC region's length) a motif-confirmed RC region is dropped
   *from the candidate pool used for the final call*, but the dropped LP
   region itself is retained separately (see step 6 below) so LucaProt's hit
   isn't lost from the record — it just doesn't win the tie-break.
   Non-overlapping LP ORFs survive independently in the normal candidate
   pool. This reuses the original overlap math (`_rc_to_nt`,
   `_overlap_fraction`), just repurposed for dropping instead of merging.
5. **Keep only the longest region per contig** — after dedup, a contig may
   still have >1 surviving candidate (e.g. two non-overlapping LP ORFs, or an
   RC region plus a non-overlapping LP ORF). Only the single longest one (by
   AA length) becomes that contig's final call. Every contig ends up tagged
   with exactly one `source`: `RdRpCATCH` or `LucaProt`.
6. **Attach LucaProt evidence, independent of who won** — after picking the
   final call, the final per-contig table gets extra columns recording
   whether LucaProt independently identified this contig at all, regardless
   of whether its region became the final call:
   - `lucaprot_also_identified` — `True` if LucaProt produced a
     motif-confirmed ORF on this contig, whether that ORF ended up being the
     final call (`source == "LucaProt"`) or was dropped by overlap dedup.
   - `lucaprot_overlap_dropped` — `True` specifically when LucaProt's hit
     was the one dropped by overlap dedup (i.e. RdRpCATCH's region won and
     became the final call instead).
   - `lucaprot_dropped_pssm_ABC`, `lucaprot_dropped_aa_length`,
     `lucaprot_dropped_protein_id`, `lucaprot_dropped_prob` — details of the
     dropped LP region (the longest one, if a contig had multiple dropped LP
     ORFs), populated only when `lucaprot_overlap_dropped == True`.

   This means a contig's identification-ability evaluation no longer depends
   on which source "won" — you can always tell whether LucaProt found it too.

### Where the final protein sequence comes from

Even though region *coordinates* come from RdRpCATCH/LucaProt, the actual
*sequence* extracted for the final FASTA comes from **palmscan's own trimmed
output** (`{sample}_rdrp_trimmed.faa` / `{sample}_lucaprot_rdrp_trimmed.faa`)
— the 150pp150-trimmed, motif-centered fragment (palmcore), not the raw
untrimmed RdRpCATCH/LucaProt protein. This was already true before the
redesign and didn't change.

## Downstream simplification

Old tier system → new 2-source design, everywhere:

| Old (8 categories/lists) | New (3 lists) |
|---|---|
| `High-confidence_RdRp`, `RDRPCatch_palmscan`, `RDRPCatch_lucaprot`, `RDRPCatch_only`, `lucaprot_palmscan`, `lucaprot_only`, `lucaprot` (union), `any_rdrp` | `RdRpCATCH`, `LucaProt`, `any_rdrp` |

Scripts rewritten to match:
- **`merge_rdrp_results.py`** — per-sample merge, motif filter + overlap dedup
  + longest-per-contig + LucaProt evidence tracking. Outputs
  `{sample}_rdrp_merged.tsv` (1 row/contig, final, with LucaProt evidence
  columns) and `{sample}_rdrp_regions.tsv` (all surviving candidates before
  longest-only selection, for inspection — no LucaProt evidence columns,
  those are contig-level only).
- **`combine_rdrp_all_samples.py`** — concatenates all samples, writes
  `RdRpCATCH.txt` / `LucaProt.txt` / `any_rdrp.txt` contig ID lists. No
  column filtering, so the new LucaProt evidence columns pass through
  untouched into `all_samples_rdrp_merged.tsv`.
- **`write_protein_id_lists.py`** — builds 2 protein-ID lists directly from
  the merged table's `source`/`aa_start`/`aa_end`/`frame` (RC) or
  `lucaprot_protein_id` (LP) columns. No more `regions_full` dependency —
  each contig has exactly one region now, so nothing else is needed.
- **`write_ictv_contig_lists.py`** — same simplification for the ICTV
  reference-sequence path (Step 9).
- **`summarize_multi_region_contigs.py`** (new) — for every contig with >1
  surviving candidate region after dedup, lists every candidate
  (source/length/strand/motif) and flags which one the longest-per-contig
  selection actually kept. Wired as rule `summarize_multi_region_contigs`
  (RNA samples, output `4_merged/multi_region_contigs_summary.tsv`) and
  `summarize_multi_region_contigs_ICTV` (ICTV reference path, output
  `ICTV/4_merged/multi_region_contigs_summary.tsv`) — same script handles
  both; it inserts a synthetic `sample="ICTV"` column when the input table
  has no `sample` column (the ICTV regions table is a single file, never
  passed through `combine_rdrp_all_samples.py`). Built to sanity-check that
  "keep longest" isn't silently biased toward one source — verified against
  real data (59 contigs with >1 candidate; longest genuinely wins regardless
  of source, no hidden RC bias).
- **`summarize_final_proteins.py`** (new) — a human-readable index of the
  two final protein FASTAs (`RdRpCATCH.faa`, `LucaProt.faa`). One row per
  contig: which FASTA it's in, the exact `protein_id` used to extract it,
  region coordinates, palmscan motif result, and the LucaProt evidence
  columns described above. Wired as rule `summarize_final_proteins` (RNA
  samples, output `8_RdRp_protein/final_proteins_summary.tsv`) and
  `summarize_final_proteins_ICTV` (ICTV path, output
  `ICTV/8_RdRp_protein/final_proteins_summary.tsv`).
- **`03_RDRP_identification.smk`** — removed `merge_high_confident` /
  `merge_high_confident_ICTV` rules, removed the 8-category
  `_ALL_PROTEIN_CATEGORIES` wildcard machinery, removed `rc_trimmed`/
  `lp_proteins` (raw, non-palmscan-verified) FASTA concatenation since
  nothing downstream needs them anymore. Only two source FASTAs remain:
  `all_samples_ps_trimmed.faa` (RdRpCATCH→palmscan) and
  `all_samples_lp_trimmed.faa` (LucaProt→palmscan). Added
  `summarize_multi_region_contigs`/`summarize_multi_region_contigs_ICTV` and
  `summarize_final_proteins`/`summarize_final_proteins_ICTV` rules, all four
  wired into `rule all` / `ictv_targets()`.

  Note: the two ICTV summary rules were initially missed when this feature
  was first added (only the RNA-sample side was wired) — added afterward
  once noticed. `merge_rdrp_results_ICTV` already called the same
  `merge_rdrp_results.py`, so `ICTV_rdrp_merged.tsv` already had the
  LucaProt evidence columns; only the summary *rules* were missing, not the
  underlying data.

## Output file map (per sample: `result/03_RDRP_identification/`)

```
4_merged/{sample}/{sample}_rdrp_merged.tsv     one row per contig (final call + LucaProt evidence)
4_merged/{sample}/{sample}_rdrp_regions.tsv    all surviving candidates (pre-longest)
4_merged/all_samples_rdrp_merged.tsv           combined, all samples
4_merged/all_samples_rdrp_regions.tsv          combined, all samples
4_merged/multi_region_contigs_summary.tsv      contigs with >1 surviving candidate, all candidates + kept flag

5_RdRp_contig/RdRpCATCH.txt / .fasta           contigs whose final region source = RdRpCATCH
5_RdRp_contig/LucaProt.txt / .fasta            contigs whose final region source = LucaProt
5_RdRp_contig/any_rdrp.txt / .fasta            union of both

8_RdRp_protein/all_samples_ps_trimmed.faa      ALL palmscan-scanned RC candidates (unfiltered pool)
8_RdRp_protein/all_samples_lp_trimmed.faa      ALL palmscan-scanned LP candidates (unfiltered pool)
8_RdRp_protein/RdRpCATCH_ids.txt / .faa        motif-confirmed, deduped, longest-per-contig (RC source)
8_RdRp_protein/LucaProt_ids.txt / .faa         motif-confirmed, deduped, longest-per-contig (LP source)
8_RdRp_protein/final_proteins_summary.tsv      lookup index for the two .faa files above, incl. LucaProt evidence
```

Same structure under `ICTV/` for the reference-sequence validation path,
including both summary files (`ICTV/4_merged/multi_region_contigs_summary.tsv`
and `ICTV/8_RdRp_protein/final_proteins_summary.tsv`).

**Sanity check identity:** `RdRpCATCH.faa` seq count + `LucaProt.faa` seq
count should always equal the number of data rows in
`all_samples_rdrp_merged.tsv` (one contig → exactly one final source →
exactly one extracted sequence). Verified: 4166 + 146 = 4312 ✓.

## `region_table` column reference

`{sample}_rdrp_regions.tsv` / `all_samples_rdrp_regions.tsv` columns
(candidate-region level, no LucaProt evidence columns here):

| column | meaning |
|---|---|
| `contig` | contig ID |
| `source` | `RdRpCATCH` or `LucaProt` |
| `pssm_ABC` | palmscan motif result for this region (`ABC`/`CAB` = confirmed) |
| `aa_start`, `aa_end` | RdRpCATCH AA region bounds (RC rows only) |
| `aa_length` | region length in AA (used to pick the longest) |
| `frame` | RdRpCATCH translation frame (RC rows only) |
| `strand` | `+`/`-`, derived from frame (RC) or ORF start/end order (LP) |
| `seq_len_aa` | full translated contig length in this frame (RC rows only; used for AA→nt conversion during overlap dedup) |
| `orf_start`, `orf_end` | LucaProt ORF nt coordinates (LP rows only) |
| `lucaprot_protein_id`, `lucaprot_prob` | LucaProt's own protein ID / confidence score (LP rows only) |

RC rows leave LP-only columns blank and vice versa — a row is never a mix of
both, since after Step 4 (dedup) each surviving region has a single,
unambiguous source.

## `all_samples_rdrp_merged.tsv` — extra LucaProt evidence columns

Only present on the **final per-contig table** (not the region table), added
in step 6 above:

| column | meaning |
|---|---|
| `lucaprot_also_identified` | `True`/`False` — did LucaProt independently produce a motif-confirmed ORF on this contig at all (won or lost the tie-break) |
| `lucaprot_overlap_dropped` | `True`/`False` — was that LP ORF specifically the one dropped by overlap dedup (RC won instead) |
| `lucaprot_dropped_pssm_ABC` | palmscan motif result of the dropped LP region, if any |
| `lucaprot_dropped_aa_length` | AA length of the dropped LP region, if any |
| `lucaprot_dropped_protein_id` | LucaProt's protein ID for the dropped region, if any |
| `lucaprot_dropped_prob` | LucaProt's confidence score for the dropped region, if any |

If a contig's final `source` is already `LucaProt`, these "dropped" columns
are all blank (nothing was dropped — LucaProt's region won outright).

## Three bugs hit and fixed during rollout

1. **Missing config keys** — `03_RDRP_identification.smk` referenced
   `rdrp_catch_output_dir`, `palmscan_output_dir`, `lucaprot_output_dir`,
   `orffinder_dir`, `rdrp_log_dir`, none of which existed in `config.yaml`.
   Added them, matching the actual existing DelftBlue directory layout
   (`orffinder_dir` ended up `0_orffinder`, not the initially-guessed
   `3_orffinder`).

2. **`EmptyDataError` on samples with zero motif-confirmed regions** —
   `merge_rdrp_results.py` called `pd.read_csv` on RdRpCATCH/palmscan/
   LucaProt per-sample files without checking if they were empty (0 bytes,
   which happens for low-diversity samples with zero hits). Added
   `os.path.getsize(...) == 0` guards before every read. Also added a fixed
   `REGION_COLUMNS` schema via `.reindex()` before writing output, so an
   empty result always writes a proper header-only TSV instead of a
   columnless blank file that crashes downstream `pd.read_csv`.

3. **`RdRpCATCH.faa` extracted 0 sequences despite 4166 IDs** —
   `write_protein_id_lists.py`'s `_rc_id()` normalized `aa_start`/`aa_end`
   via `int(float(...))` (correctly handling pandas' `27` → `"27.0"` float
   coercion after `pd.concat`/`reindex` mixed int and NaN columns), but
   **did not** apply the same normalization to `frame` — so the built ID
   contained `frame=-1.0` while the actual FASTA header said `frame=-1`,
   and seqkit grep matched nothing. Fixed by applying `int(float(...))`
   to `frame` too.

## Empirical check: is "keep longest" biased toward one source?

Ran `summarize_multi_region_contigs` against real DelftBlue data: only ~59
contigs (out of 4312) had more than one surviving candidate region after
dedup. Spot-checking several: the longest region wins regardless of source
— RC wins when it's genuinely longer, LP wins when *it's* genuinely longer
(e.g. `AnQing_U_32_0000000065`: RC=309aa lost to LP=355aa). No hidden RC
bias in the final selection — the only place RdRpCATCH gets an automatic
edge is the overlap-dedup tie-break (step 4), not the longest-per-contig
pick (step 5).

## Step 10: nr viral-origin check (Diamond BLASTp), added 2026-07-14

Motivation (from the RVMT/Neri-style methods text the user is matching):
after motif-based filtering, candidates can still include endogenous viral
elements (EVEs) or false positives whose closest relatives in nr are
cellular/host proteins rather than RNA virus proteins. A final Diamond
BLASTp-against-nr check, keeping only sequences with a significant RNA
virus hit, removes these before the candidates are called final RdRps.

**Design choice:** rather than resolving each hit's taxon lineage in Python
(originally implemented with `ete3.NCBITaxa`, then dropped per user
feedback — one more conda env/dependency to maintain for something diamond
already does natively), the taxon check is pushed entirely into diamond
itself via `--taxonlist <riboviria_taxid>` at search time. The Diamond
database is built with `--taxonmap` + `--taxonnodes` + `--taxonnames` (NCBI
accession2taxid + taxdump), which lets `diamond blastp --taxonlist 2559587`
restrict the search to hits whose taxon falls anywhere under **Riboviria**
(NCBI taxid `2559587`, the realm containing all RNA viruses) — diamond
walks the taxonomic tree internally. A query passes the filter if and only
if that restricted search returns a hit within the e-value threshold; no
separate taxonomy library is needed anywhere downstream.

### What was added

- **`config/config.yaml`** — new `diamond_nr:` block: `db_dir`, `nr_fasta`,
  `dmnd_db`, `taxdump_dir` (+ `taxdump_url`, needed at *build* time for
  `--taxonnodes`/`--taxonnames`), `prot_accession2taxid`, `db_marker`,
  `riboviria_taxid: 2559587` (passed to `diamond blastp --taxonlist` at
  search time), `evalue` (default `0.001`, matching NCBI nr BLASTp
  convention), plus standard `threads`/`runtime`/`memory`/`partition`/
  `account` keys (`memory` partition, since nr is large). **Paths are
  placeholders** — point `db_dir`/`dmnd_db`/`taxdump_dir` at the real
  DelftBlue locations once built (or let `00_setup_RNA_database.smk`
  download+build them from scratch).

  **NCBI FASTA deprecation note:** NCBI stopped publishing the plain
  `nr.gz` FASTA on the FTP site in April 2024
  ([announcement](https://ncbiinsights.ncbi.nlm.nih.gov/2024/01/25/blast-fasta-unavailable-on-ftp/)).
  Only the pre-formatted BLAST database is still served, so `nr_fasta` and
  `prot_accession2taxid` are **regenerated locally** from that BLAST db via
  `blastdbcmd`, not downloaded directly (see below) — `nr_fasta` is now a
  plain (non-gzipped) FASTA path.
- **`workflow/00_setup_RNA_database.smk`** — two new rules (internet-
  dependent downloads only; run on the login node, matching this file's
  existing "no sbatch, no compute-node internet" convention):
  - `download_blast_nr_db` — runs `update_blastdb.pl --decompress nr`
    (ships with `blast+`, same tool NCBI recommends post-deprecation) to
    fetch and assemble the multi-volume pre-formatted BLAST protein
    database into `diamond_nr.db_dir`. Writes a `.blastdb_download.done`
    marker (the BLAST db is many volume files, not one clean output path).
  - `download_diamond_nr_taxdump` — downloads + extracts NCBI's
    `new_taxdump.tar.gz` (`nodes.dmp`/`names.dmp`, needed to build the
    taxonomic tree diamond uses for `--taxonlist`).
  - Both wired into `00_setup_RNA_database.smk`'s `rule all`.
- **`workflow/00b_build_diamond_nr_db.smk`** (new file) — two rules that
  never touch the network, so they're split out to run via `sbatch` on a
  compute node instead of the login node, after the two download rules
  above have finished:
  - `generate_nr_fasta_and_taxid_map` — regenerates both the FASTA
    (`blastdbcmd -db nr -entry all -out nr_fasta`) and an
    accession→taxid mapping (`blastdbcmd -db nr -entry all -outfmt "%a %T"`)
    from the already-downloaded BLAST db in one pass each — no separate
    `prot.accession2taxid.FULL.gz` download needed, since the BLAST db
    already carries taxid info per record. `blastdbcmd` streams
    sequentially (doesn't load nr into RAM) and is single-threaded; the
    real cost here is disk (FASTA + taxid map add ~100+ GB alongside the
    BLAST db already on disk) and wall time, not memory.
  - `build_diamond_nr_db` — `diamond makedb --taxonmap <acc2tax>
    --taxonnodes <nodes.dmp> --taxonnames <names.dmp>`, so `diamond blastp`
    can both report per-hit `staxids` and filter searches by
    `--taxonlist <taxid>` directly. Writes a `db_marker` touch-file (same
    convention as `checkv`/`genomad`/`virsorter2` db rules) since the
    actual `.dmnd` file is huge and not a good Snakemake output target by
    itself.
  - Run it with e.g.
    `snakemake --snakefile workflow/00b_build_diamond_nr_db.smk --use-conda
    --cores 16 --rerun-triggers input` submitted via `sbatch`.
- **`scripts/03_RDRP_identification/filter_nr_viral_origin.py`** (new,
  no ete3/taxonomy-library dependency — uses only pandas) — reads the
  diamond output, which the calling rule already restricted to
  `--taxonlist <riboviria_taxid> --evalue <threshold> -k 1`, i.e. it only
  ever contains, per query, that query's single best hit *within Riboviria*
  (if one exists at all). A query with no row in that file has no
  significant RNA virus match in nr and is dropped. Writes a filtered
  FASTA (`{pcat}_viral_confirmed.faa`) and a full audit table
  (`{pcat}_nr_filter_table.tsv`, one row per query: kept/dropped, hit
  accession, staxids, e-value, reason).
- **`workflow/03_RDRP_identification.smk`** — new Step 10, applied to both
  final protein FASTAs (`RdRpCATCH.faa`, `LucaProt.faa`) via a
  `{pcat}` wildcard (`pcat = "RdRpCATCH|LucaProt"`):
  - `rule diamond_blastp_nr_riboviria` — `diamond blastp -k 1
    --taxonlist {config[diamond_nr][riboviria_taxid]} --evalue
    {config[diamond_nr][evalue]} --outfmt 6 qseqid sseqid evalue bitscore
    staxids` against the nr `.dmnd` db.
  - `rule filter_nr_viral_origin` — runs the script above (env:
    `../envs/python.yaml`, no new env needed), output under
    `result/03_RDRP_identification/9_nr_viral_check/`.
  - Mirrored for the ICTV reference path as
    `diamond_blastp_nr_riboviria_ICTV` / `filter_nr_viral_origin_ICTV`,
    output under `ICTV/9_nr_viral_check/` (`_ICTV_NR_DIR`, defined
    alongside the other `_ICTV_*_DIR` path constants near the top of the
    file).
  - Both wired into `rule all` (RNA path) and `ictv_targets()` (ICTV path).

### Output files (new)

```
9_nr_viral_check/{pcat}_diamond_nr_riboviria.tsv   diamond hits restricted to Riboviria (one row per query with a hit)
9_nr_viral_check/{pcat}_viral_confirmed.faa        final filtered FASTA -- only sequences with a Riboviria hit survive
9_nr_viral_check/{pcat}_nr_filter_table.tsv        audit table: every query, kept/dropped, hit id, taxids, reason
```

Same structure under `ICTV/9_nr_viral_check/` for the reference path.
`{pcat}` is `RdRpCATCH` or `LucaProt`, matching the Step 8 FASTA names.

### Still to do before running on DelftBlue

The `diamond_nr` config block currently has placeholder paths
(`database/diamond_nr/...`). Either:
- run `workflow/00_setup_RNA_database.smk` (login node) to download the
  BLAST nr db + taxdump, then `workflow/00b_build_diamond_nr_db.smk` (via
  `sbatch`, compute node) to regenerate the FASTA/taxid map and build the
  Diamond db (large — 100+ GB, hours to build, needs the `memory`
  partition), or
- if an nr Diamond database already exists somewhere on DelftBlue, point
  `diamond_nr.dmnd_db` at it directly instead of re-downloading — but it
  must have been built with `--taxonmap` **and** `--taxonnodes`/
  `--taxonnames` (not just `--taxonmap` alone) for `--taxonlist` filtering
  to work at search time; otherwise rebuild it with `build_diamond_nr_db`
  in `workflow/00b_build_diamond_nr_db.smk`.
