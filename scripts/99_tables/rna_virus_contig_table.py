"""
Build RNA virus contig abundance table (3_RNA_virus_table).

Steps:
  1. Load contig RPKMF and read_count tables; keep only unknown::acc: rows.
     Strip prefix to get bare contig_id. Extract Asm_length from read_count.
  2. Join RdRp info from all_samples_rdrp_merged.tsv on contig_id.
     If multiple RdRp rows per contig (same contig, different ORF), keep the
     row with the longest ORF (aa_length).
     source values: "RdRpCATCH" or "LucaProt".
  3. Join bait info from bait_target_map.tsv on contig_id.
     If multiple bait rows per contig, keep highest confidence:
     RdRpCATCH > LucaProt > RVMT.
  4. Add contig_tag column combining RdRp + bait evidence:
     RdRpCATCH > LucaProt > bait_RdRpCATCH > bait_LucaProt > bait_RVMT > unknown.
  5. For contigs still tagged "unknown", look up GeNomad's
     Unified Genome Composition (08_other_virus_identification annotated summary)
     by contig_id == seq_name. Rows with composition "ssRNA" are skipped (already
     covered by RdRp). Non-ssRNA hits get tag "genomad_<composition>"
     (e.g. "genomad_ssDNA", "genomad_dsDNA"). Never overrides an existing tag.
  6. Reorder: info cols | Asm_length | RPKMF sample cols | read_count sample cols.

Usage:
    python scripts/99_tables/rna_virus_contig_table.py
"""

import os
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────
RPKMF_TSV        = "result/99_tables/01_cp_table/04_modified_esvirtue_ribodetector/esvirtu_rpkmf.subspecies.tsv"
COUNT_TSV        = "result/99_tables/01_cp_table/04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv"
RDRP_MERGED_TSV  = "result/99_tables/01_cp_table/03_RDRP_identification/all_samples_rdrp_merged.tsv"
BAIT_TSV         = "result/99_tables/01_cp_table/06_contig_bait/bait_target_map.tsv"
GENOMAD_TSV      = "result/99_tables/01_cp_table/08_other_virus_identification/all_samples_virus_filter_summary.annotated.tsv"
OUT_TSV          = "result/99_tables/3_RNA_virus_table/rna_virus_contig_table.tsv"

# ── source priority (lower index = higher priority) ───────────────────────────
SOURCE_ORDER = ["RdRpCATCH", "LucaProt"]
BAIT_ORDER   = ["RdRpCATCH", "LucaProt", "RVMT"]
TAG_ORDER    = ["RdRpCATCH", "LucaProt",
                "bait_RdRpCATCH", "bait_LucaProt", "bait_RVMT",
                "unknown"]

# composition values already covered by RdRp — skip in GeNomad fallback
GENOMAD_SKIP_COMPOSITIONS = {"ssRNA"}


def source_rank(source, order):
    try:
        return order.index(source)
    except ValueError:
        return len(order)


# =============================================================================

def load_contigs(path):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "subspecies"})
    before = len(df)
    df = df[df["subspecies"].str.startswith("unknown::acc:")].reset_index(drop=True)
    print(f"  Kept {len(df)} contig rows, removed {before - len(df)} ESvirtu reference rows.")
    df["contig_id"] = df["subspecies"].str.removeprefix("unknown::acc:")
    df = df.drop(columns=["subspecies"])
    cols = ["contig_id"] + [c for c in df.columns if c != "contig_id"]
    return df[cols]


def load_count_contigs(path):
    """Load read_count table; return (count_df, asm_length_series) for contig rows only."""
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "subspecies"})
    df = df[df["subspecies"].str.startswith("unknown::acc:")].reset_index(drop=True)
    df["contig_id"] = df["subspecies"].str.removeprefix("unknown::acc:")
    asm_length = df.set_index("contig_id")["Asm_length"]
    count_cols = [c for c in df.columns if c not in ("subspecies", "contig_id", "Asm_length")]
    count_df = df[["contig_id"] + count_cols]
    return count_df, asm_length


def build_rdrp_lookup(path):
    """Best RdRp hit per contig from all_samples_rdrp_merged.tsv.

    Tie-breaking: if both RdRpCATCH and LucaProt detected the same contig,
    keep RdRpCATCH (higher source_rank priority); within the same source,
    keep the row with the longest ORF (aa_length).
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["aa_length"] = pd.to_numeric(df.get("aa_length"), errors="coerce").fillna(0)
    df["_src_rank"] = df["source"].apply(lambda s: source_rank(s, SOURCE_ORDER))
    df = (
        df.sort_values(["_src_rank", "aa_length"], ascending=[True, False])
          .drop_duplicates("contig", keep="first")
    )
    keep_cols = ["contig", "source", "pssm_ABC", "aa_length", "strand"]
    keep_cols = [c for c in keep_cols if c in df.columns]
    return df[keep_cols].set_index("contig")


def build_bait_lookup(path):
    """Best-confidence bait row per contig, including source info."""
    bait = pd.read_csv(path, sep="\t")
    bait["_rank"] = bait["bait_category"].apply(lambda t: source_rank(t, BAIT_ORDER))
    bait = bait.sort_values("_rank").drop_duplicates("recovered_contig", keep="first")
    return bait[["recovered_contig", "bait_category", "best_bait_contig", "bait_source",
                 "pident", "tcov", "evalue", "n_bait_hits"]].set_index("recovered_contig")


def build_genomad_lookup(path):
    """contig_id (seq_name) -> Unified Genome Composition, excluding ssRNA rows."""
    gm = pd.read_csv(path, sep="\t", usecols=["seq_name", "Unified Genome Composition"])
    gm = gm.rename(columns={"seq_name": "contig_id"})
    gm = gm[~gm["Unified Genome Composition"].isin(GENOMAD_SKIP_COMPOSITIONS)]
    gm = gm[gm["Unified Genome Composition"].notna() & (gm["Unified Genome Composition"] != "")]
    gm = gm.drop_duplicates("contig_id", keep="first")
    return gm.set_index("contig_id")["Unified Genome Composition"]


def assign_tag(row):
    rdrp_source = row.get("source", pd.NA)
    bait_cat    = row.get("bait_category", pd.NA)

    # Own RdRp contig (RdRpCATCH or LucaProt) takes top priority
    if pd.notna(rdrp_source) and rdrp_source in ("RdRpCATCH", "LucaProt"):
        return rdrp_source

    # Bait-recovered contig
    if pd.notna(bait_cat) and str(bait_cat).strip():
        return f"bait_{bait_cat}"

    return "unknown"


# =============================================================================

def main():
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

    print("Loading contig RPKMF table...")
    df = load_contigs(RPKMF_TSV)
    rpkmf_sample_cols = [c for c in df.columns if c != "contig_id"]

    print("Loading contig read_count table...")
    count_df, asm_length = load_count_contigs(COUNT_TSV)
    count_sample_cols = [c for c in count_df.columns if c != "contig_id"]

    print("Adding Asm_length...")
    df = df.join(asm_length, on="contig_id", how="left")

    print("Building RdRp lookup...")
    rdrp = build_rdrp_lookup(RDRP_MERGED_TSV)
    df = df.join(rdrp, on="contig_id", how="left")

    print("Building bait lookup...")
    bait = build_bait_lookup(BAIT_TSV)
    df = df.join(bait, on="contig_id", how="left")

    print("Assigning contig_tag...")
    df["contig_tag"] = df.apply(assign_tag, axis=1)

    print("Applying GeNomad tags to remaining unknown contigs...")
    genomad = build_genomad_lookup(GENOMAD_TSV)
    unknown_mask = df["contig_tag"] == "unknown"
    genomad_hits = df.loc[unknown_mask, "contig_id"].map(genomad)
    genomad_tag = "genomad_" + genomad_hits
    df.loc[unknown_mask, "contig_tag"] = genomad_tag.combine_first(df.loc[unknown_mask, "contig_tag"])
    print(f"  Tagged {genomad_hits.notna().sum()} previously-unknown contigs from GeNomad.")

    print("Joining read_count sample columns...")
    count_idx = count_df.set_index("contig_id")[count_sample_cols]
    df = df.join(count_idx, on="contig_id", how="left")

    # ── column order ──────────────────────────────────────────────────────────
    info_cols = ["contig_id", "Asm_length", "source", "pssm_ABC", "aa_length", "strand",
                 "bait_category", "best_bait_contig", "bait_source",
                 "pident", "tcov", "evalue", "n_bait_hits", "contig_tag"]
    info_cols = [c for c in info_cols if c in df.columns]
    df = df[info_cols + rpkmf_sample_cols + count_sample_cols]

    # ── summary ───────────────────────────────────────────────────────────────
    counts = df["contig_tag"].value_counts()
    known = counts.reindex(TAG_ORDER, fill_value=0)
    genomad_tags = counts[~counts.index.isin(TAG_ORDER)].sort_index()
    print("\ncontig_tag summary:")
    print(known.to_string())
    if not genomad_tags.empty:
        print(genomad_tags.to_string())
    print()

    print(f"Writing {OUT_TSV} ({len(df)} rows, {len(df.columns)} columns)...")
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print("Done.")


if __name__ == "__main__":
    main()
