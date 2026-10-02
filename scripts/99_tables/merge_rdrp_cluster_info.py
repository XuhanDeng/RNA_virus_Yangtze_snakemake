"""
Merge protein/contig/cluster-level seqkit fx2tab outputs into one summary
table under result/99_tables/3_RDRP_Summary/.

Inputs (all from result/99_tables/01_cp_table/03_RDRP_identification/):
  10_final/nr_filtered/combined_full.fx2tab.tsv       -- per-protein: seq_id, length, GC
  10_final/nr_filtered/rdrp_contigs.fx2tab.tsv        -- per-contig:  seq_id, length, GC
  10_final/nr_filtered_cluster/combined_full_c90.faa.clstr -- CD-HIT cluster membership
  21_taxonomic_assignment_cluster/full_length/full_length_diamond_rvmt_tophit.tsv
      -- DIAMOND-vs-RVMT taxonomy, searched only against the CD-HIT cluster
      representatives (Step 21 runs on the same combined_full_c90.faa set),
      so qseqid here corresponds to representative_seq_id, not every protein.
  10_final/final_rdrp_merged.tsv
      -- final NR-confirmed per-contig RdRp identification table (source
      tool, motif region coordinates, LucaProt cross-evidence); one row per
      contig_id, joined directly (not via representative_seq_id -- every
      contig has its own independent call here, regardless of clustering).

Output (result/99_tables/3_RDRP_Summary/):
  nr_filtered_protein_contig_cluster.tsv
      One row per NR-filtered protein (all ~4153 of them, not just cluster
      representatives):
        protein_length                   -- from the protein itself
        contig_id, contig_length, contig_GC -- its parent contig (joined via
            the same protein-ID -> contig-ID transform used elsewhere in
            this pipeline: strip ORF\\d+_ prefix, truncate at first ':',
            strip trailing _frame=...)
        cluster_id                       -- CD-HIT cluster this protein fell
            into at Step 20's 90% AAI clustering
        representative_seq_id            -- the seq_id of the '*'-marked
            representative for that protein's cluster
        is_cluster_representative        -- True if this row's own seq_id IS
            that representative (audit column; True for exactly one row per
            cluster_id)
        cluster_size                     -- how many total sequences fell
            into this protein's cluster
        rvmt_* columns (sseqid, pident, length, evalue, bitscore, qcovhsp,
            scovhsp, hit_rank, Phylum, Class, Order, Family, Genus,
            confidence_tier) -- the representative's own DIAMOND-vs-RVMT
            taxonomy call, joined on representative_seq_id (NOT this row's
            own seq_id) so every protein in a cluster shares its
            representative's annotation.
        rdrp_* columns (source, pssm_ABC, aa_start, aa_end, aa_length, frame,
            strand, seq_len_aa, orf_start, orf_end, lucaprot_protein_id,
            lucaprot_prob, lucaprot_also_identified, lucaprot_overlap_dropped,
            lucaprot_dropped_pssm_ABC, lucaprot_dropped_aa_length,
            lucaprot_dropped_protein_id, lucaprot_dropped_prob) -- from
            final_rdrp_merged.tsv, joined on contig_id (this row's own
            contig, independent of clustering).

Usage:
    python scripts/99_tables/merge_rdrp_cluster_info.py
"""

import os
import re
import pandas as pd

# ── paths ─────────────────────────────────────────────────────────────────────
_BASE      = "result/99_tables/01_cp_table/03_RDRP_identification/10_final"
_TAX_BASE  = "result/99_tables/01_cp_table/03_RDRP_identification/21_taxonomic_assignment_cluster"

PROTEIN_TSV   = f"{_BASE}/nr_filtered/combined_full.fx2tab.tsv"
CONTIG_TSV    = f"{_BASE}/nr_filtered/rdrp_contigs.fx2tab.tsv"
CLSTR_FILE    = f"{_BASE}/nr_filtered_cluster/combined_full_c90.faa.clstr"
TOPHIT_TSV    = f"{_TAX_BASE}/full_length/full_length_diamond_rvmt_tophit.tsv"
RDRP_MERGED_TSV = f"{_BASE}/final_rdrp_merged.tsv"

# tophit.tsv columns to carry over, renamed with an rvmt_ prefix (qseqid is
# consumed as the join key, not carried over as its own column).
_TOPHIT_COLS = ["sseqid", "pident", "length", "evalue", "bitscore",
                "qcovhsp", "scovhsp", "hit_rank", "Phylum", "Class",
                "Order", "Family", "Genus", "confidence_tier"]

# final_rdrp_merged.tsv columns to carry over, renamed with an rdrp_ prefix
# (contig is consumed as the join key, not carried over as its own column).
_RDRP_MERGED_COLS = [
    "source", "pssm_ABC", "aa_start", "aa_end", "aa_length", "frame",
    "strand", "seq_len_aa", "orf_start", "orf_end", "lucaprot_protein_id",
    "lucaprot_prob", "lucaprot_also_identified", "lucaprot_overlap_dropped",
    "lucaprot_dropped_pssm_ABC", "lucaprot_dropped_aa_length",
    "lucaprot_dropped_protein_id", "lucaprot_dropped_prob",
]

OUT_DIR            = "result/99_tables/3_RDRP_Summary"
OUT_PROTEIN_CONTIG = f"{OUT_DIR}/nr_filtered_protein_contig_cluster.tsv"


def to_contig(pid: str) -> str:
    """Same transform used in make_final_output.py: protein/ORF ID -> parent contig ID."""
    pid = re.sub(r"^ORF\d+_", "", pid)
    pid = pid.split(":")[0]
    pid = re.sub(r"_frame=.*", "", pid)
    return pid


def parse_clstr(path):
    """
    Parse a CD-HIT .clstr file.

    Returns:
      cluster_of[seq_id]         -> cluster_id (int)
      representative_of[seq_id]  -> True/False (is this seq_id the '*'-marked
                                     representative of its cluster?)
      members_of[cluster_id]     -> list of every seq_id in that cluster, in
                                     .clstr file order
      rep_seq_of[cluster_id]     -> seq_id of that cluster's representative
    """
    cluster_of = {}
    representative_of = {}
    members_of = {}
    rep_seq_of = {}
    current_cluster = None

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">Cluster"):
                current_cluster = int(line.split()[-1])
                members_of[current_cluster] = []
                continue
            # e.g. "0\t1514aa, >ORF4_WuHu_D_37_0000000010:4547:3... *"
            #      "5\t888aa, >WuHan_D_26_0000000113_frame=3_RdRp_146-1033... at 100.00%"
            m = re.search(r">(\S+?)\.\.\.\s*(\*|at\s)", line)
            if not m:
                continue
            seq_id = m.group(1)
            is_rep = m.group(2).strip() == "*"
            cluster_of[seq_id] = current_cluster
            representative_of[seq_id] = is_rep
            members_of[current_cluster].append(seq_id)
            if is_rep:
                rep_seq_of[current_cluster] = seq_id

    return cluster_of, representative_of, members_of, rep_seq_of


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    print("Loading protein and contig fx2tab tables...")
    protein_df = pd.read_csv(PROTEIN_TSV, sep="\t")
    contig_df  = pd.read_csv(CONTIG_TSV, sep="\t")

    print("Loading RVMT tophit taxonomy table...")
    tophit_df = pd.read_csv(TOPHIT_TSV, sep="\t")
    tophit_df = tophit_df[["qseqid"] + _TOPHIT_COLS]
    tophit_df = tophit_df.rename(columns={"qseqid": "representative_seq_id"})
    tophit_df = tophit_df.rename(columns={c: f"rvmt_{c}" for c in _TOPHIT_COLS})
    print(f"  {len(tophit_df)} representative(s) with RVMT taxonomy.")

    print("Loading final_rdrp_merged.tsv (per-contig RdRp identification)...")
    rdrp_merged_df = pd.read_csv(RDRP_MERGED_TSV, sep="\t")
    rdrp_merged_df = rdrp_merged_df[["contig"] + _RDRP_MERGED_COLS]
    rdrp_merged_df = rdrp_merged_df.rename(columns={"contig": "contig_id"})
    rdrp_merged_df = rdrp_merged_df.rename(columns={c: f"rdrp_{c}" for c in _RDRP_MERGED_COLS})
    print(f"  {len(rdrp_merged_df)} contig(s) with RdRp identification info.")

    print("Parsing CD-HIT .clstr file...")
    cluster_of, representative_of, members_of, rep_seq_of = parse_clstr(CLSTR_FILE)
    print(f"  Parsed {len(cluster_of)} sequences across "
          f"{len(members_of)} clusters.")

    print("Building protein/contig/cluster merge...")
    protein_df = protein_df.rename(columns={"length": "protein_length"}).drop(columns=["GC"])
    protein_df["contig_id"] = protein_df["seq_id"].apply(to_contig)

    contig_df = contig_df.rename(columns={"seq_id": "contig_id",
                                            "length": "contig_length",
                                            "GC": "contig_GC"})

    merged = protein_df.merge(contig_df, on="contig_id", how="left")
    merged["cluster_id"] = merged["seq_id"].map(cluster_of)
    merged["representative_seq_id"] = merged["cluster_id"].map(rep_seq_of)
    merged["is_cluster_representative"] = merged["seq_id"].map(representative_of)
    merged["cluster_size"] = merged["cluster_id"].map(
        lambda cid: len(members_of[cid]) if pd.notna(cid) and cid in members_of else pd.NA
    )

    print("Joining RVMT taxonomy by representative_seq_id (same annotation "
          "shared across every protein in a cluster)...")
    merged = merged.merge(tophit_df, on="representative_seq_id", how="left")

    print("Joining final_rdrp_merged.tsv by contig_id (independent per contig, "
          "not shared across a cluster)...")
    merged = merged.merge(rdrp_merged_df, on="contig_id", how="left")

    n_missing_contig = merged["contig_length"].isna().sum()
    n_missing_cluster = merged["cluster_id"].isna().sum()
    n_missing_tax = merged["rvmt_Phylum"].isna().sum()
    n_missing_rdrp = merged["rdrp_source"].isna().sum()
    if n_missing_contig:
        print(f"  WARNING: {n_missing_contig} protein(s) had no matching contig_id "
              f"in {CONTIG_TSV}.")
    if n_missing_cluster:
        print(f"  WARNING: {n_missing_cluster} protein(s) had no cluster_id "
              f"in {CLSTR_FILE}.")
    if n_missing_tax:
        print(f"  WARNING: {n_missing_tax} protein(s) had no RVMT taxonomy "
              f"for their cluster's representative in {TOPHIT_TSV}.")
    if n_missing_rdrp:
        print(f"  WARNING: {n_missing_rdrp} protein(s) had no matching row "
              f"in {RDRP_MERGED_TSV}.")

    n_reps = int(merged["is_cluster_representative"].sum())
    n_clusters = len(members_of)
    if n_reps != n_clusters:
        print(f"  WARNING: found {n_reps} representative rows but "
              f"{n_clusters} clusters -- expected exactly one representative "
              f"per cluster.")
    else:
        print(f"  Confirmed exactly one representative per cluster "
              f"({n_reps} representatives across {n_clusters} clusters).")

    print(f"  Cluster sizes: min={merged['cluster_size'].min()}, "
          f"max={merged['cluster_size'].max()}, "
          f"mean={merged['cluster_size'].mean():.2f}")

    rvmt_cols = [f"rvmt_{c}" for c in _TOPHIT_COLS]
    rdrp_cols = [f"rdrp_{c}" for c in _RDRP_MERGED_COLS]
    merged = merged[["seq_id", "contig_id", "protein_length",
                      "contig_length", "contig_GC", "cluster_id",
                      "representative_seq_id", "is_cluster_representative",
                      "cluster_size"] + rvmt_cols + rdrp_cols]
    merged.to_csv(OUT_PROTEIN_CONTIG, sep="\t", index=False)
    print(f"Wrote {OUT_PROTEIN_CONTIG} ({len(merged)} rows, "
          f"covering all NR-filtered proteins).")

    print("Done.")


if __name__ == "__main__":
    main()
