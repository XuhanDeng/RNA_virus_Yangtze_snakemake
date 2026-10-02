"""
Recalculate RPKMF and TPM for RNA viruses only
(5_Recalculated_RPKM_result/0_RDRP_Esvirtue).

Combines two sources:
  - esvirtu_info_table.tsv   -- ESvirtu-detected known viruses; kept rows are
                                 those with genome_type in {ssRNA(+), ssRNA(-), dsRNA}.
  - rna_virus_contig_table.tsv -- de novo contigs; kept rows are those tagged
                                 RdRpCATCH or LucaProt (genomad_* and unknown
                                 rows are excluded -- not RNA virus evidence).

Steps:
  1. Load esvirtu_info_table.tsv; keep rows where genome_type in
     {ssRNA(+), ssRNA(-), dsRNA}. Drop old *_RPKMF columns.
  2. Load rna_virus_contig_table.tsv; keep rows where contig_tag is
     RdRpCATCH or LucaProt.
  3. Rename row-ID columns to seq_id; add source column; concat into one table.
     All metadata columns from both tables are retained (NA where absent).
  4. Recalculate RPKMF per sample using the combined RNA-virus read counts:
       RPKMF = (count / Asm_length / 1000) / (total_counts_in_sample / 1e6)
  5. Recalculate TPM per sample (two-pass normalization):
       rpk = count / (Asm_length / 1000)
       TPM = rpk / (sum(rpk in sample) / 1e6)
  6. Write three output files:
       rna_virus_recalc_read_count.tsv  — combined table with raw read counts
       rna_virus_recalc_rpkmf.tsv       — same but sample cols replaced by RPKMF
       rna_virus_recalc_tpm.tsv         — same but sample cols replaced by TPM

RVMT taxonomy is NOT joined here -- see annotate_rna_virus_taxonomy.py, which
reads these three outputs and writes taxonomy-annotated copies to
5_Recalculated_RPKM_result/1_taxonmy/.

Usage:
    python scripts/99_tables/rna_virus_recalc_rpkmf.py
"""

import os
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────
ESVIRTU_TSV   = "result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv"
CONTIG_TSV    = "result/99_tables/4_contig_tag_abundance/rna_virus_contig_table.tsv"
OUTDIR        = "result/99_tables/5_Recalculated_RPKM_result/0_RDRP_Esvirtue"
OUT_COUNT     = OUTDIR + "/rna_virus_recalc_read_count.tsv"
OUT_RPKMF     = OUTDIR + "/rna_virus_recalc_rpkmf.tsv"
OUT_TPM       = OUTDIR + "/rna_virus_recalc_tpm.tsv"

# ── filters ───────────────────────────────────────────────────────────────────
RNA_GENOME_TYPES = {"ssRNA(+)", "ssRNA(-)", "dsRNA"}
RNA_CONTIG_TAGS  = {"RdRpCATCH", "LucaProt"}


# =============================================================================

def load_esvirtu(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    before = len(df)
    df = df[df["genome_type"].isin(RNA_GENOME_TYPES)].reset_index(drop=True)
    print(f"  ESvirtu: kept {len(df)} RNA-virus rows (removed {before - len(df)}).")

    # drop old RPKMF columns — we will recalculate
    rpkmf_cols = [c for c in df.columns if c.endswith("_RPKMF")]
    df = df.drop(columns=rpkmf_cols)

    df = df.rename(columns={"subspecies": "seq_id"})
    df.insert(1, "source", "esvirtu")
    return df


def load_contigs(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    before = len(df)
    df = df[df["contig_tag"].isin(RNA_CONTIG_TAGS)].reset_index(drop=True)
    print(f"  Contigs: kept {len(df)} RNA-virus rows (removed {before - len(df)}).")

    df = df.rename(columns={"contig_id": "seq_id"})
    if "source" not in df.columns:
        df.insert(1, "source", "contig")
    return df


def recalc_rpkmf(combined, count_cols):
    """Return a copy of combined with count_cols replaced by recalculated RPKMF."""
    counts = combined[count_cols].fillna(0).values.astype(float)
    lengths = combined["Asm_length"].fillna(1).values.astype(float)

    # reads per kilobase: count / (length / 1000)
    rpk = counts / (lengths[:, None] / 1000.0)

    # per-sample total mapped to RNA virus set (in millions)
    totals = counts.sum(axis=0) / 1e6   # shape (n_samples,)
    totals[totals == 0] = 1             # avoid div-by-zero for empty samples

    rpkmf_values = rpk / totals[None, :]

    rpkmf_cols = [c.replace("_read_count", "_rpkmf") for c in count_cols]
    rpkmf = pd.DataFrame(rpkmf_values, columns=rpkmf_cols, index=combined.index)

    result = combined.drop(columns=count_cols).copy()
    result = pd.concat([result, rpkmf], axis=1)
    return result, rpkmf_cols


def recalc_tpm(combined, count_cols):
    """Return a copy of combined with count_cols replaced by recalculated TPM."""
    counts = combined[count_cols].fillna(0).values.astype(float)
    lengths = combined["Asm_length"].fillna(1).values.astype(float)

    # reads per kilobase: count / (length / 1000)
    rpk = counts / (lengths[:, None] / 1000.0)

    # per-sample sum of RPK (in millions) — the TPM normalization denominator
    rpk_totals = rpk.sum(axis=0) / 1e6   # shape (n_samples,)
    rpk_totals[rpk_totals == 0] = 1      # avoid div-by-zero for empty samples

    tpm_values = rpk / rpk_totals[None, :]

    tpm_cols = [c.replace("_read_count", "_tpm") for c in count_cols]
    tpm = pd.DataFrame(tpm_values, columns=tpm_cols, index=combined.index)

    result = combined.drop(columns=count_cols).copy()
    result = pd.concat([result, tpm], axis=1)
    return result, tpm_cols


# =============================================================================

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("Loading ESvirtu RNA virus rows...")
    esvirtu = load_esvirtu(ESVIRTU_TSV)

    print("Loading RNA virus contig rows (RdRpCATCH, LucaProt only)...")
    contigs = load_contigs(CONTIG_TSV)

    print("Combining tables...")
    combined = pd.concat([esvirtu, contigs], axis=0, ignore_index=True)
    print(f"  Combined: {len(combined)} rows total.")

    # identify sample read_count columns (same set in both tables)
    count_cols = [c for c in combined.columns if c.endswith("_read_count")]
    print(f"  Sample columns: {len(count_cols)}")

    # ── column order: seq_id | source | all metadata | Asm_length | sample cols
    meta_cols = [c for c in combined.columns
                 if c not in {"seq_id", "source", "Asm_length"} and c not in count_cols]
    col_order = ["seq_id", "source"] + meta_cols + ["Asm_length"] + count_cols
    col_order = [c for c in col_order if c in combined.columns]
    combined = combined[col_order]

    # ── QC: per-sample total read counts ─────────────────────────────────────
    print("\nPer-sample total RNA virus read counts:")
    totals = combined[count_cols].fillna(0).sum(axis=0)
    for col, val in totals.items():
        print(f"  {col}: {int(val):,}")
    print()

    # ── write read_count output ───────────────────────────────────────────────
    print(f"Writing {OUT_COUNT} ({len(combined)} rows, {len(combined.columns)} cols)...")
    combined.to_csv(OUT_COUNT, sep="\t", index=False)

    # ── recalculate RPKMF ─────────────────────────────────────────────────────
    print("Recalculating RPKMF...")
    rpkmf_df, rpkmf_cols = recalc_rpkmf(combined, count_cols)

    # reorder: same meta + Asm_length, then rpkmf cols
    non_count = [c for c in col_order if c not in count_cols]
    rpkmf_df = rpkmf_df[non_count + rpkmf_cols]

    print(f"Writing {OUT_RPKMF} ({len(rpkmf_df)} rows, {len(rpkmf_df.columns)} cols)...")
    rpkmf_df.to_csv(OUT_RPKMF, sep="\t", index=False)

    # ── recalculate TPM ────────────────────────────────────────────────────────
    print("Recalculating TPM...")
    tpm_df, tpm_cols = recalc_tpm(combined, count_cols)
    tpm_df = tpm_df[non_count + tpm_cols]

    print(f"Writing {OUT_TPM} ({len(tpm_df)} rows, {len(tpm_df.columns)} cols)...")
    tpm_df.to_csv(OUT_TPM, sep="\t", index=False)

    print("Done.")


if __name__ == "__main__":
    main()
