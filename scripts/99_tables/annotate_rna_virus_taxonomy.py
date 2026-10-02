"""
Annotate the recalculated RNA-virus abundance tables with taxonomy
(5_Recalculated_RPKM_result/1_taxonmy).

Reads the three tables produced by rna_virus_recalc_rpkmf.py
(0_RDRP_Esvirtue/) and:
  1. Joins rvmt_* RVMT DIAMOND taxonomy columns from
     nr_filtered_protein_contig_cluster.tsv onto contig rows (seq_id == its
     contig_id there). ESvirtu rows already carry their own taxonomy (Name,
     kingdom, phylum, tclass, order, family, genus, species, subspecies,
     TaxID, etc. from esvirtu_info_table.tsv) -- those columns are left
     untouched; esvirtu rows simply get NA for the rvmt_* columns, since RVMT
     taxonomy is a contig-only annotation never computed for ESvirtu
     reference hits.
  2. Adds five unified tax_* columns (tax_phylum, tax_class, tax_order,
     tax_family, tax_genus) so both sources can be read from one place
     regardless of which pipeline produced them:
       - esvirtu rows -> from phylum / tclass / order / family / genus,
                         with the GTDB-style rank prefix (p__/c__/o__/f__/g__)
                         stripped so values match RVMT's plain-name format
                         (e.g. "p__Duplornaviricota" -> "Duplornaviricota").
       - contig rows  -> from rvmt_Phylum / rvmt_Class / rvmt_Order /
                         rvmt_Family / rvmt_Genus (already unprefixed).
     The original esvirtu_* and rvmt_* columns are kept as-is, unprefixed
     values included -- only the new unified tax_* columns are cleaned.

Steps:
  1. Load nr_filtered_protein_contig_cluster.tsv; keep contig_id + rvmt_*
     columns only. One row per contig is assumed (no duplicate RdRp proteins
     per contig); a warning is printed if that assumption doesn't hold.
  2. For each of read_count / rpkmf / tpm, join rvmt_* onto seq_id == contig_id,
     then fill tax_* from whichever source (esvirtu or rvmt) is populated.
  3. Write annotated copies to 1_taxonmy/.

Usage:
    python scripts/99_tables/annotate_rna_virus_taxonomy.py
"""

import os
import re
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────
RECALC_DIR   = "result/99_tables/5_Recalculated_RPKM_result/0_RDRP_Esvirtue"
IN_COUNT     = RECALC_DIR + "/rna_virus_recalc_read_count.tsv"
IN_RPKMF     = RECALC_DIR + "/rna_virus_recalc_rpkmf.tsv"
IN_TPM       = RECALC_DIR + "/rna_virus_recalc_tpm.tsv"

TAXONOMY_TSV = "result/99_tables/3_RDRP_Summary/nr_filtered_protein_contig_cluster.tsv"

OUTDIR       = "result/99_tables/5_Recalculated_RPKM_result/1_taxonmy"
OUT_COUNT    = OUTDIR + "/rna_virus_recalc_read_count.annotated.tsv"
OUT_RPKMF    = OUTDIR + "/rna_virus_recalc_rpkmf.annotated.tsv"
OUT_TPM      = OUTDIR + "/rna_virus_recalc_tpm.annotated.tsv"

# rvmt_* taxonomy columns to join in from nr_filtered_protein_contig_cluster.tsv
_RVMT_COLS = ["rvmt_sseqid", "rvmt_pident", "rvmt_length", "rvmt_evalue",
              "rvmt_bitscore", "rvmt_qcovhsp", "rvmt_scovhsp", "rvmt_hit_rank",
              "rvmt_Phylum", "rvmt_Class", "rvmt_Order", "rvmt_Family",
              "rvmt_Genus", "rvmt_confidence_tier"]

# unified tax_* column -> (esvirtu source column, rvmt source column)
_TAX_UNIFY = {
    "tax_phylum": ("phylum", "rvmt_Phylum"),
    "tax_class":  ("tclass", "rvmt_Class"),
    "tax_order":  ("order",  "rvmt_Order"),
    "tax_family": ("family", "rvmt_Family"),
    "tax_genus":  ("genus",  "rvmt_Genus"),
}


# =============================================================================

def load_taxonomy(path):
    """contig_id -> rvmt_* taxonomy columns from nr_filtered_protein_contig_cluster.tsv.

    One row per contig is assumed (no duplicate RdRp proteins per contig);
    a warning is printed if that assumption doesn't hold.
    """
    df = pd.read_csv(path, sep="\t", low_memory=False)
    keep_cols = ["contig_id"] + [c for c in _RVMT_COLS if c in df.columns]
    df = df[keep_cols]
    n_before = len(df)
    df = df.drop_duplicates("contig_id", keep="first")
    if len(df) != n_before:
        print(f"  WARNING: {n_before - len(df)} duplicate contig_id row(s) in "
              f"{path}, kept first occurrence of each.")
    return df.set_index("contig_id")


def _strip_rank_prefix(series):
    """Strip a leading GTDB-style rank prefix (p__, c__, o__, f__, g__, ...)."""
    return series.astype("string").str.replace(r"^[a-z]__", "", regex=True)


def add_unified_tax_columns(df):
    """Add tax_phylum/tax_class/tax_order/tax_family/tax_genus, filled from
    whichever source (esvirtu or rvmt) is populated for each row. ESvirtu's
    values have their GTDB-style rank prefix stripped so they match RVMT's
    plain-name format."""
    for tax_col, (esvirtu_col, rvmt_col) in _TAX_UNIFY.items():
        esvirtu_series = df[esvirtu_col] if esvirtu_col in df.columns else pd.Series(pd.NA, index=df.index)
        esvirtu_series = _strip_rank_prefix(esvirtu_series)
        rvmt_series    = df[rvmt_col]    if rvmt_col    in df.columns else pd.Series(pd.NA, index=df.index)
        df[tax_col] = esvirtu_series.combine_first(rvmt_series)
    return df


def annotate(in_path, out_path, taxonomy):
    df = pd.read_csv(in_path, sep="\t", low_memory=False)
    df = df.join(taxonomy, on="seq_id", how="left")
    n_tax = df["rvmt_Phylum"].notna().sum() if "rvmt_Phylum" in df.columns else 0
    print(f"  {os.path.basename(in_path)}: matched RVMT taxonomy for "
          f"{n_tax}/{len(df)} rows.")

    df = add_unified_tax_columns(df)
    n_unified = df["tax_phylum"].notna().sum()
    print(f"  {os.path.basename(in_path)}: tax_phylum populated for "
          f"{n_unified}/{len(df)} rows.")

    print(f"Writing {out_path} ({len(df)} rows, {len(df.columns)} cols)...")
    df.to_csv(out_path, sep="\t", index=False)


# =============================================================================

def main():
    os.makedirs(OUTDIR, exist_ok=True)

    print("Loading RVMT taxonomy...")
    taxonomy = load_taxonomy(TAXONOMY_TSV)

    print("Annotating recalculated RNA-virus tables...")
    annotate(IN_COUNT, OUT_COUNT, taxonomy)
    annotate(IN_RPKMF, OUT_RPKMF, taxonomy)
    annotate(IN_TPM, OUT_TPM, taxonomy)

    print("Done.")


if __name__ == "__main__":
    main()
