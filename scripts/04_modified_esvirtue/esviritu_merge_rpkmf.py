#!/usr/bin/env python3
"""
Merge per-sample EsViritu assembly_summary TSVs and compute RPKMF wide matrices.

Replaces two former scripts (esviritu_rpkmf.py + esviritu_merge_rpkmf.py).

Input:  per-sample .merged.tsv files (raw, with read_count / Asm_length /
        filtered_reads_in_sample columns).

Outputs (derived from --output path, e.g. .../all_samples...rpkmf.tsv):
  assembly level:
    <output>                          wide RPKMF matrix  (assembly rows)
    <base>.read_count.tsv             wide read_count matrix with metadata

  subspecies / species level:
    <base>.subspecies.tsv             wide RPKMF  (sum reads → recompute RPKMF)
    <base>.species.tsv                wide RPKMF
    <base_plain>.subspecies.read_count.tsv
    <base_plain>.species.read_count.tsv

RPKMF formula (all levels):
    RPKMF = (read_count / (total_Asm_length_kb)) / (filtered_reads / 1e6)

For subspecies/species the denominator uses the SUM of Asm_length across all
contigs that belong to the same label (per sample), so the unit stays RPKMF.
"""

import argparse
import pandas as pd


def _pick_column(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


def _compute_rpkmf(reads: pd.Series, length: pd.Series, filtered_reads: float) -> pd.Series:
    """RPKMF = reads / (length_kb) / (filtered_reads_per_million)."""
    if filtered_reads == 0:
        return pd.Series(0.0, index=reads.index)
    rpkmf = (reads / (length / 1000.0)) / (filtered_reads / 1e6)
    return rpkmf.where(length > 0, 0.0).fillna(0.0)


def _is_unknown_label(series: pd.Series) -> pd.Series:
    s = series.fillna("").astype(str).str.strip()
    lc = s.str.lower()
    return (
        s.eq("")
        | lc.isin(["unknow", "unknown", "nan", "none", "na"])
        | lc.str.startswith("unknown::")
        | lc.str.startswith("unknow::")
    )


def main():
    parser = argparse.ArgumentParser(
        description="Merge per-sample EsViritu TSVs and compute RPKMF wide matrices."
    )
    parser.add_argument("--inputs", nargs="+", required=True,
                        help="Per-sample .merged.tsv files (raw, with read_count)")
    parser.add_argument("--output", required=True,
                        help="Assembly-level wide RPKMF TSV path")
    args = parser.parse_args()

    # ------------------------------------------------------------------ #
    # Pass 1: load each per-sample TSV, compute RPKMF, collect metadata   #
    # ------------------------------------------------------------------ #
    per_sample = {}     # sample_name -> DataFrame (all original columns + _rpkmf)
    meta_cols_union = set()
    id_col_name = None

    for path in args.inputs:
        df = pd.read_csv(path, sep="\t", dtype=str)
        # drop repeated header rows from cat-based merges
        if len(df.columns) > 0:
            df = df[df[df.columns[0]] != df.columns[0]]

        id_col = _pick_column(df.columns, ["Assembly", "Accession"])
        if id_col is None:
            raise SystemExit(f"No Assembly/Accession column in {path}")
        if id_col_name is None:
            id_col_name = id_col
        elif id_col_name != id_col:
            df = df.rename(columns={id_col: id_col_name})

        sample_col = _pick_column(df.columns, ["sample_ID", "sample", "Sample"])
        if sample_col is None:
            raise SystemExit(f"No sample_ID column in {path}")
        sample_vals = df[sample_col].dropna().unique()
        sample = str(sample_vals[0]) if len(sample_vals) > 0 else "sample"

        filtered_col = _pick_column(df.columns, ["filtered_reads_in_sample", "filtered_reads"])
        length_col   = _pick_column(df.columns, ["Asm_length", "Asm_Length", "Assembly_length",
                                                  "assembly_length", "Length", "Contig_length"])
        read_col     = _pick_column(df.columns, ["read_count", "Read_count", "reads"])

        missing = [n for n, v in [("filtered_reads_in_sample", filtered_col),
                                   ("Asm_length", length_col),
                                   ("read_count", read_col)] if v is None]
        if missing:
            raise SystemExit(f"{path}: missing required columns: {', '.join(missing)}")

        per_sample_cols = {
            sample_col,
            filtered_col,
            "covered_bases",
            "avg_read_identity",
        }
        meta_cols = [c for c in df.columns if c not in per_sample_cols]
        meta_cols_union.update(meta_cols)

        # compute per-row RPKMF
        filtered_vals = pd.to_numeric(df[filtered_col].dropna(), errors="coerce").dropna()
        filtered_reads = float(filtered_vals.iloc[0]) if not filtered_vals.empty else 0.0

        length = pd.to_numeric(df[length_col], errors="coerce").fillna(0.0)
        reads  = pd.to_numeric(df[read_col],   errors="coerce").fillna(0.0)

        df[length_col] = length
        df[read_col]   = reads
        df["_rpkmf"]   = _compute_rpkmf(reads, length, filtered_reads)
        df["_filtered_reads"] = filtered_reads

        per_sample[sample] = df

    # ------------------------------------------------------------------ #
    # Build metadata table: one row per assembly, first-non-empty wins    #
    # ------------------------------------------------------------------ #
    meta_parts = []
    for df in per_sample.values():
        keep = [c for c in df.columns if c in meta_cols_union]
        meta_parts.append(df[keep])
    meta_all = pd.concat(meta_parts, ignore_index=True)
    for c in meta_cols_union:
        if c not in meta_all.columns:
            meta_all[c] = pd.NA

    def _first_nonempty(series):
        for v in series:
            if pd.notna(v) and str(v).strip() != "":
                return v
        return pd.NA

    agg_cols = {c: _first_nonempty for c in meta_cols_union if c != id_col_name}
    meta_df = meta_all.groupby(id_col_name, as_index=False).agg(agg_cols)

    # relabel unknown subspecies/species with unique assembly-based labels
    assembly_col = "Assembly" if "Assembly" in meta_df.columns else id_col_name
    for tax_col in ["subspecies", "species"]:
        if tax_col not in meta_df.columns:
            continue
        is_unk = _is_unknown_label(meta_df[tax_col])
        meta_df.loc[is_unk, tax_col] = "unknown::" + meta_df.loc[is_unk, assembly_col].astype(str)

    # ------------------------------------------------------------------ #
    # Assembly-level wide matrices                                         #
    # ------------------------------------------------------------------ #
    rpkmf_wide = None
    rc_wide    = None

    for sample, df in per_sample.items():
        read_col = _pick_column(df.columns, ["read_count", "Read_count", "reads"])

        per = df[[id_col_name, "_rpkmf"]].copy().rename(columns={"_rpkmf": f"{sample}_RPKMF"})
        rpkmf_wide = per if rpkmf_wide is None else pd.merge(rpkmf_wide, per, on=id_col_name, how="outer")

        per_rc = df[[id_col_name, read_col]].copy().rename(columns={read_col: f"{sample}_read_count"})
        rc_wide = per_rc if rc_wide is None else pd.merge(rc_wide, per_rc, on=id_col_name, how="outer")

    rpkmf_wide = rpkmf_wide.fillna(0)
    rc_wide    = rc_wide.fillna(0)

    merged_rpkmf = pd.merge(meta_df, rpkmf_wide, on=id_col_name, how="outer")
    merged_rc    = pd.merge(meta_df, rc_wide,    on=id_col_name, how="outer")

    # output path helpers
    base       = args.output[:-4] if args.output.endswith(".tsv") else args.output
    base_plain = base[:-6]        if base.endswith(".rpkmf")      else base

    merged_rpkmf.to_csv(args.output, sep="\t", index=False)
    merged_rc.to_csv(base_plain + ".read_count.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ #
    # Aggregation: sum reads + length per taxon label, then recompute     #
    # RPKMF — ensures counts are summed before normalisation              #
    # ------------------------------------------------------------------ #
    length_col_meta = _pick_column(list(meta_df.columns),
                                   ["Asm_length", "Asm_Length", "Assembly_length",
                                    "assembly_length", "Length", "Contig_length"])

    def _aggregate_level(tax_col, out_rpkmf_path, out_rc_path):
        if tax_col not in merged_rc.columns:
            return

        rc_sample_cols = [c for c in merged_rc.columns if c.endswith("_read_count")]
        if not rc_sample_cols:
            return

        cols_needed = [tax_col, assembly_col] + rc_sample_cols
        has_length = bool(length_col_meta and length_col_meta in merged_rc.columns)
        if has_length:
            merged_rc[length_col_meta] = pd.to_numeric(
                merged_rc[length_col_meta], errors="coerce"
            ).fillna(0.0)
            cols_needed.append(length_col_meta)

        df_work = merged_rc[cols_needed].copy()
        for c in rc_sample_cols:
            df_work[c] = pd.to_numeric(df_work[c], errors="coerce").fillna(0.0)

        # unknown rows get a unique per-assembly label
        is_unk = _is_unknown_label(df_work[tax_col])
        df_work.loc[is_unk, tax_col] = "unknown::" + df_work.loc[is_unk, assembly_col].astype(str)

        # sum read_counts per label; use longest contig length as the RPKMF denominator
        rc_summed = df_work.groupby(tax_col, as_index=False)[rc_sample_cols].sum()
        if has_length:
            max_length = df_work.groupby(tax_col, as_index=False)[length_col_meta].max()
            rc_summed = pd.merge(rc_summed, max_length, on=tax_col, how="left")
        rc_summed.to_csv(out_rc_path, sep="\t", index=False)

        # recompute RPKMF: summed reads / longest contig length
        rpkmf_rows = {tax_col: rc_summed[tax_col].values}
        for sample, df in per_sample.items():
            rc_col_name = f"{sample}_read_count"
            if rc_col_name not in rc_summed.columns:
                continue
            fr    = float(df["_filtered_reads"].iloc[0])
            reads = rc_summed[rc_col_name]
            if has_length and length_col_meta in rc_summed.columns:
                length = rc_summed[length_col_meta]
            else:
                # no length info: fall back to RPM
                length = pd.Series(1000.0, index=reads.index)
            rpkmf_rows[f"{sample}_RPKMF"] = _compute_rpkmf(reads, length, fr).values

        pd.DataFrame(rpkmf_rows).to_csv(out_rpkmf_path, sep="\t", index=False)

    _aggregate_level("subspecies", base + ".subspecies.tsv",
                     base_plain + ".subspecies.read_count.tsv")
    _aggregate_level("species",    base + ".species.tsv",
                     base_plain + ".species.read_count.tsv")


if __name__ == "__main__":
    main()
