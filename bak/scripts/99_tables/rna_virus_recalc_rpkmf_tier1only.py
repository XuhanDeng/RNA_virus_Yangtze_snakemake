"""
Recalculate RPKMF for RNA viruses — tier1 + bait_tier1 only.

Same logic as rna_virus_recalc_rpkmf.py but keeps only tier1a, tier1b,
bait_tier1a, bait_tier1b (excludes tier2, tier3, bait_tier2, bait_RVMT).
Outputs go to 5_Recalculated_RPKM_result/3_RNA_virus_Recalculated_RPKM_tier1only/.
"""

import os
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────
ESVIRTU_TSV   = "result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv"
CONTIG_TSV    = "result/99_tables/3_RNA_virus_table/rna_virus_contig_table.tsv"
OUTDIR        = "result/99_tables/5_Recalculated_RPKM_result/3_RNA_virus_Recalculated_RPKM_tier1only"
OUT_COUNT     = OUTDIR + "/rna_virus_recalc_read_count.tsv"
OUT_RPKMF     = OUTDIR + "/rna_virus_recalc_rpkmf.tsv"
OUT_TPM       = OUTDIR + "/rna_virus_recalc_tpm.tsv"

# ── filters (tier1 + bait_tier1 only) ────────────────────────────────────────
RNA_GENOME_TYPES = {"ssRNA(+)", "ssRNA(-)", "dsRNA"}

RNA_CONTIG_TAGS = {
    "tier1a", "tier1b",
    "bait_tier1a", "bait_tier1b",
}


# =============================================================================

def load_esvirtu(path):
    df = pd.read_csv(path, sep="\t", low_memory=False)
    before = len(df)
    df = df[df["genome_type"].isin(RNA_GENOME_TYPES)].reset_index(drop=True)
    print(f"  ESvirtu: kept {len(df)} RNA-virus rows (removed {before - len(df)}).")

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

    rpkmf_cols = [c for c in df.columns if c.endswith("_RPKMF")]
    df = df.drop(columns=rpkmf_cols)

    df = df.rename(columns={"contig_id": "seq_id"})
    df.insert(1, "source", "contig")
    return df


def recalc_rpkmf(combined, count_cols):
    counts = combined[count_cols].fillna(0).values.astype(float)
    lengths = combined["Asm_length"].fillna(1).values.astype(float)

    rpk = counts / (lengths[:, None] / 1000.0)

    totals = counts.sum(axis=0) / 1e6
    totals[totals == 0] = 1

    rpkmf_values = rpk / totals[None, :]

    rpkmf_cols = [c.replace("_read_count", "_rpkmf") for c in count_cols]
    rpkmf = pd.DataFrame(rpkmf_values, columns=rpkmf_cols, index=combined.index)

    result = combined.drop(columns=count_cols).copy()
    result = pd.concat([result, rpkmf], axis=1)
    return result, rpkmf_cols


def recalc_tpm(combined, count_cols):
    counts = combined[count_cols].fillna(0).values.astype(float)
    lengths = combined["Asm_length"].fillna(1).values.astype(float)

    rpk = counts / (lengths[:, None] / 1000.0)

    rpk_totals = rpk.sum(axis=0) / 1e6
    rpk_totals[rpk_totals == 0] = 1

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

    print("Loading RNA virus contig rows (tier1 only: tier1a/1b, bait_tier1a/1b)...")
    contigs = load_contigs(CONTIG_TSV)

    print("Combining tables...")
    combined = pd.concat([esvirtu, contigs], axis=0, ignore_index=True)
    print(f"  Combined: {len(combined)} rows total.")

    count_cols = [c for c in combined.columns if c.endswith("_read_count")]
    print(f"  Sample columns: {len(count_cols)}")

    meta_cols = [c for c in combined.columns
                 if c not in {"seq_id", "source", "Asm_length"} and c not in count_cols]
    col_order = ["seq_id", "source"] + meta_cols + ["Asm_length"] + count_cols
    col_order = [c for c in col_order if c in combined.columns]
    combined = combined[col_order]

    print("\nPer-sample total RNA virus read counts (tier1 only):")
    totals = combined[count_cols].fillna(0).sum(axis=0)
    for col, val in totals.items():
        print(f"  {col}: {int(val):,}")
    print()

    print(f"Writing {OUT_COUNT} ({len(combined)} rows, {len(combined.columns)} cols)...")
    combined.to_csv(OUT_COUNT, sep="\t", index=False)

    print("Recalculating RPKMF...")
    rpkmf_df, rpkmf_cols = recalc_rpkmf(combined, count_cols)

    non_count = [c for c in col_order if c not in count_cols]
    rpkmf_df = rpkmf_df[non_count + rpkmf_cols]

    print(f"Writing {OUT_RPKMF} ({len(rpkmf_df)} rows, {len(rpkmf_df.columns)} cols)...")
    rpkmf_df.to_csv(OUT_RPKMF, sep="\t", index=False)

    print("Recalculating TPM...")
    tpm_df, tpm_cols = recalc_tpm(combined, count_cols)
    tpm_df = tpm_df[non_count + tpm_cols]

    print(f"Writing {OUT_TPM} ({len(tpm_df)} rows, {len(tpm_df.columns)} cols)...")
    tpm_df.to_csv(OUT_TPM, sep="\t", index=False)

    print("Done.")


if __name__ == "__main__":
    main()
