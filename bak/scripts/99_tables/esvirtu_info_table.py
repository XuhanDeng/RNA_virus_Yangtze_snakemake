"""
Build enriched ESvirtu abundance table (2_esvirtu_table).

Steps:
  1. Filter RPKMF and read_count tables to ESviritu rows only (t__ prefix).
     Extract Asm_length from read_count table (last column).
  2. Join with virus_pathogen_database.all_metadata.tsv on subspecies →
     taxonomy columns + Accession.
  3. Add genome_type column from phylum mapping.
  4. Join on TaxID with virushostdb.tsv → vhd host columns
     (exclude entries where host tax id = 1 / "root").
  5. Join on Accession with host_metadata.tsv → host_meta_Host column.
  6. Add merged_host: vhd_host_name first; fall back to host_meta_Host.
  7. Add minimum_host_lineage: for viruses with multiple hosts, the lowest
     common ancestor (LCA) across all host lineages in vhd_host_lineage —
     the most specific taxon shared by every host of that virus.
  8. Reorder: info cols | Asm_length | RPKMF sample cols | read_count sample cols.
  9. Print genome_type summary counts.

Usage:
    python scripts/99_tables/esvirtu_info_table.py
"""

import os
import warnings
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────
RPKMF_TSV       = "result/99_tables/01_cp_table/04_modified_esvirtue_ribodetector/esvirtu_rpkmf.subspecies.tsv"
COUNT_TSV       = "result/99_tables/01_cp_table/04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv"
ALL_META_TSV    = "database/esviritu/v3.2.4/virus_pathogen_database.all_metadata.tsv"
VIRUSHOSTDB_TSV = "database/virushostdb/virushostdb.tsv"
HOST_META_TSV   = "database/esviritu/v3.2.4/virus_pathogen_database.host_metadata.v3.2.4.tsv"
OUT_TSV         = "result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv"

# ── taxonomy columns to pull from all_metadata ────────────────────────────────
META_COLS = ["Name", "Segment", "kingdom", "phylum", "tclass", "order",
             "family", "genus", "species", "subspecies", "TaxID"]

# ── virushostdb columns to add ────────────────────────────────────────────────
VHD_COLS = ["virus tax id", "virus name", "host tax id", "host name",
            "host lineage", "pmid", "evidence"]

# ── phylum → genome_type mapping (complete ICTV 2023 phyla) ──────────────────
PHYLUM_TO_GENOME_TYPE = {
    # ssDNA
    "p__Cressdnaviricota":      "ssDNA",   # circoviruses, nanoviruses
    "p__Cossaviricota":         "ssDNA",   # parvoviruses, densoviruses
    "p__Preplasmiviricota":     "ssDNA",   # anelloviruses, genomoviruses, smacoviruses
    # ssRNA(+)
    "p__Pisuviricota":          "ssRNA(+)",  # picornaviruses, coronaviruses, potyviruses
    "p__Kitrinoviricota":       "ssRNA(+)",  # flaviviruses, tombusviruses, hypoviruses
    "p__Lenarviricota":         "ssRNA(+)",  # mitoviruses, narnaviruses, leviviruses
    # ssRNA(-)
    "p__Negarnaviricota":       "ssRNA(-)",  # rhabdo, bunya, paramyxo, orthomyxo
    # dsRNA
    "p__Duplornaviricota":      "dsRNA",   # reoviruses, partitiviruses, totiviruses
    # dsDNA
    "p__Nucleocytoviricota":    "dsDNA",   # giant viruses (mimiviruses, poxviruses)
    "p__Peploviricota":         "dsDNA",   # herpesviruses
    "p__Uroviricota":           "dsDNA",   # tailed bacteriophages
    "p__Dividoviricota":        "dsDNA",   # tectiviruses, corticoviruses
    "p__Helvetiaviricota":      "dsDNA",   # salterprovirus
    # dsDNA-RT
    "p__Artverviricota":        "dsDNA-RT",  # retroviruses, hepadnaviruses
    # ssDNA-RT (caulimoviruses)
    "p__Shotokuvirae":          "ssDNA-RT",  # note: Shotokuvirae is a kingdom; phylum = Cressdnaviricota
    # unclassified
    "p__unclassified_Viruses":  "unknown",
}


# =============================================================================

def assign_genome_type(phylum_series):
    return phylum_series.map(PHYLUM_TO_GENOME_TYPE).fillna("unknown")


def load_rpkmf(path):
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "subspecies"})
    before = len(df)
    df = df[df["subspecies"].str.startswith("t__")].reset_index(drop=True)
    print(f"  Kept {len(df)} ESviritu rows, removed {before - len(df)} assembly contigs.")
    return df


def load_count(path):
    """Load read_count table; return (count_df_no_asm_length, asm_length_series)."""
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={df.columns[0]: "subspecies"})
    df = df[df["subspecies"].str.startswith("t__")].reset_index(drop=True)
    asm_length = df.set_index("subspecies")["Asm_length"]
    count_cols = [c for c in df.columns if c not in ("subspecies", "Asm_length")]
    count_df = df[["subspecies"] + count_cols]
    return count_df, asm_length


def build_metadata_lookup(path):
    """Return per-subspecies aggregated metadata from all_metadata.tsv."""
    meta = pd.read_csv(path, sep="\t", low_memory=False)

    acc_map = (meta.groupby("subspecies")["Accession"]
                   .apply(lambda x: ";".join(x.dropna().astype(str).unique()))
                   .rename("Accession"))

    rows = []
    for subsp, grp in meta.groupby("subspecies"):
        row = {"subspecies": subsp}
        for col in META_COLS:
            if col == "subspecies":
                continue
            vals = grp[col].dropna().unique()
            if len(vals) > 1:
                warnings.warn(
                    f"[WARNING] subspecies '{subsp}' has multiple distinct "
                    f"values for '{col}': {vals.tolist()} — taking first."
                )
            row[col] = vals[0] if len(vals) > 0 else pd.NA
        rows.append(row)

    tax_df = pd.DataFrame(rows).set_index("subspecies")
    tax_df = tax_df.join(acc_map)
    return tax_df


def build_virushostdb_lookup(path):
    """Return per-TaxID aggregated host info; exclude root (host tax id == 1)."""
    vhd = pd.read_csv(path, sep="\t", low_memory=False)
    vhd = vhd.rename(columns={"virus tax id": "virus_tax_id"})

    # drop rows where host is root (tax id == 1 or host name == "root")
    if "host tax id" in vhd.columns:
        vhd = vhd[vhd["host tax id"].astype(str) != "1"]
    if "host name" in vhd.columns:
        vhd = vhd[vhd["host name"].str.strip().str.lower() != "root"]

    rows = []
    for taxid, grp in vhd.groupby("virus_tax_id"):
        row = {"TaxID": taxid}
        for col in VHD_COLS:
            if col == "virus tax id":
                continue
            if col not in vhd.columns:
                continue
            vals = grp[col].dropna().unique()
            row[f"vhd_{col.replace(' ', '_')}"] = ";".join(vals.astype(str)) if len(vals) else pd.NA
        rows.append(row)

    return pd.DataFrame(rows).set_index("TaxID")


def build_host_meta_lookup(path):
    """Return per-Accession Host from host_metadata tsv; exclude root entries."""
    hm = pd.read_csv(path, sep="\t", low_memory=False)
    if "Host" not in hm.columns:
        return pd.Series(dtype=str, name="host_meta_Host")
    # drop rows where Host is "root", "NA", or empty
    hm = hm[~hm["Host"].str.strip().str.lower().isin(["root", "na", ""])]
    return (hm.groupby("Accession")["Host"]
              .apply(lambda x: ";".join(x.dropna().astype(str).unique()))
              .rename("host_meta_Host"))


def merge_host(row):
    """Use vhd_host_name if available, else host_meta_Host."""
    vhd = row.get("vhd_host_name", pd.NA)
    meta = row.get("host_meta_Host", pd.NA)
    if pd.notna(vhd) and str(vhd).strip():
        return vhd
    return meta if pd.notna(meta) else pd.NA


LINEAGE_ROOT = "cellular organisms"


def split_host_lineages(lineage_field: str) -> list:
    """
    vhd_host_lineage concatenates one full lineage per host with no separator
    between them — each lineage simply starts again with "cellular organisms".
    Split back into a list of per-host lineage lists (each a list of taxa).
    """
    text = str(lineage_field)
    # split before each occurrence of the root taxon (except a leading one)
    parts = text.split(LINEAGE_ROOT + ";")
    lineages = []
    for i, part in enumerate(parts):
        if not part.strip() and i == 0:
            continue
        taxa = [LINEAGE_ROOT] + [t.strip() for t in part.split(";") if t.strip()]
        lineages.append(taxa)
    return lineages


def lowest_common_ancestor(lineage_field) -> str:
    """Full lineage path (root -> most specific shared taxon) common to every host."""
    if pd.isna(lineage_field) or not str(lineage_field).strip():
        return pd.NA
    lineages = split_host_lineages(lineage_field)
    if not lineages:
        return pd.NA
    if len(lineages) == 1:
        return "; ".join(lineages[0])
    common = []
    for taxon_a, *rest in zip(*lineages):
        if all(taxon_a == taxon_b for taxon_b in rest):
            common.append(taxon_a)
        else:
            break
    return "; ".join(common) if common else pd.NA


# =============================================================================

def main():
    os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

    print("Loading RPKMF table...")
    rpkmf = load_rpkmf(RPKMF_TSV)
    rpkmf_sample_cols = [c for c in rpkmf.columns if c != "subspecies"]

    print("Loading read_count table...")
    count_df, asm_length = load_count(COUNT_TSV)
    count_sample_cols = [c for c in count_df.columns if c != "subspecies"]

    print("Building metadata lookup...")
    tax_df = build_metadata_lookup(ALL_META_TSV)

    print("Joining taxonomy...")
    merged = rpkmf.join(tax_df, on="subspecies", how="left")

    print("Adding Asm_length...")
    merged = merged.join(asm_length, on="subspecies", how="left")

    print("Adding genome_type...")
    merged["genome_type"] = assign_genome_type(merged["phylum"])

    print("Building virushostdb lookup...")
    vhd_df = build_virushostdb_lookup(VIRUSHOSTDB_TSV)
    merged["TaxID"] = pd.to_numeric(merged["TaxID"], errors="coerce")
    merged = merged.join(vhd_df, on="TaxID", how="left")

    print("Building host metadata lookup...")
    host_series = build_host_meta_lookup(HOST_META_TSV)
    merged["_acc_first"] = merged["Accession"].str.split(";").str[0].str.strip()
    merged = merged.join(host_series, on="_acc_first", how="left")
    merged = merged.drop(columns=["_acc_first"])

    print("Merging host columns (vhd priority)...")
    merged["merged_host"] = merged.apply(merge_host, axis=1)

    print("Computing minimum_host_lineage (LCA across hosts)...")
    merged["minimum_host_lineage"] = merged["vhd_host_lineage"].apply(lowest_common_ancestor)

    print("Joining read_count sample columns...")
    count_idx = count_df.set_index("subspecies")[count_sample_cols]
    merged = merged.join(count_idx, on="subspecies", how="left")

    # ── column order ──────────────────────────────────────────────────────────
    vhd_out_cols = [f"vhd_{c.replace(' ', '_')}" for c in VHD_COLS if c != "virus tax id"]
    info_cols = (["subspecies"] +
                 [c for c in META_COLS if c != "subspecies"] +
                 ["genome_type", "Asm_length", "Accession"] +
                 vhd_out_cols +
                 ["host_meta_Host", "merged_host", "minimum_host_lineage"])
    info_cols = [c for c in info_cols if c in merged.columns]
    merged = merged[info_cols + rpkmf_sample_cols + count_sample_cols]

    # ── genome_type summary ───────────────────────────────────────────────────
    print("\ngenome_type summary:")
    print(merged["genome_type"].value_counts().to_string())
    print()

    print(f"Writing {OUT_TSV} ({len(merged)} rows, {len(merged.columns)} columns)...")
    merged.to_csv(OUT_TSV, sep="\t", index=False)
    print("Done.")


if __name__ == "__main__":
    main()
