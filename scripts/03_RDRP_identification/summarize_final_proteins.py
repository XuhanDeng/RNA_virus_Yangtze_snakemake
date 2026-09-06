#!/usr/bin/env python3
"""
Build a human-readable index of the two final protein FASTAs
(RdRpCATCH.faa, LucaProt.faa) in 8_RdRp_protein/.

For every contig in the final merged table, show which source it came from
(RdRpCATCH or LucaProt), the exact protein ID used to extract its sequence
from the corresponding FASTA, its region coordinates, its palmscan motif
result (pssm_ABC), and whether LucaProt ALSO independently identified this
contig -- even when RdRpCATCH's region won the overlap tie-break and became
the final call. Useful for evaluating LucaProt's identification ability
independent of which region source ended up "winning" per contig.

Reads:
  --merged  all_samples_rdrp_merged.tsv (or ICTV_rdrp_merged.tsv)

Writes:
  --output  one row per contig, with a `protein_id` column matching the
            header in RdRpCATCH.faa / LucaProt.faa (minus the ">").
"""

import argparse
import re
import sys

import pandas as pd


def _lp_label(protein_id: str) -> str:
    s = str(protein_id).strip().lstrip(">")
    s = re.sub(r"^lcl\|", "", s)
    return s.split()[0]


def _rc_id(row) -> str:
    try:
        frame = str(int(float(row.get("frame", ""))))
        aa_s  = str(int(float(row.get("aa_start", ""))))
        aa_e  = str(int(float(row.get("aa_end", ""))))
    except (TypeError, ValueError):
        return None
    return f"{row['contig']}_frame={frame}_RdRp_{aa_s}-{aa_e}"


def _lp_id(row) -> str:
    pid = row.get("lucaprot_protein_id", "")
    if pd.notna(pid) and str(pid).strip() not in ("", "nan"):
        return _lp_label(str(pid))
    return None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--merged", required=True,
                        help="all_samples_rdrp_merged.tsv (or ICTV_rdrp_merged.tsv)")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    df = pd.read_csv(args.merged, sep="\t", dtype=str)
    print(f"Loaded {len(df)} contigs from {args.merged}", file=sys.stderr)

    def _protein_id(row):
        source = str(row.get("source", "")).strip()
        if source == "RdRpCATCH":
            return _rc_id(row)
        if source == "LucaProt":
            return _lp_id(row)
        return None

    df["protein_id"] = df.apply(_protein_id, axis=1)
    df["fasta_file"] = df["source"].map({"RdRpCATCH": "RdRpCATCH.faa", "LucaProt": "LucaProt.faa"})

    cols = [c for c in ["sample", "contig", "source", "fasta_file", "protein_id",
                        "pssm_ABC", "aa_length", "strand",
                        "aa_start", "aa_end", "frame",
                        "orf_start", "orf_end",
                        "lucaprot_protein_id", "lucaprot_prob",
                        "lucaprot_also_identified", "lucaprot_overlap_dropped",
                        "lucaprot_dropped_pssm_ABC", "lucaprot_dropped_aa_length",
                        "lucaprot_dropped_protein_id", "lucaprot_dropped_prob"]
            if c in df.columns]
    out = df[cols]
    out.to_csv(args.output, sep="\t", index=False)

    print(f"Wrote {len(out)} rows -> {args.output}", file=sys.stderr)
    print(f"  by source: {out['source'].value_counts().to_dict()}", file=sys.stderr)
    n_missing = out["protein_id"].isna().sum()
    if n_missing:
        print(f"  WARNING: {n_missing} rows have no resolvable protein_id", file=sys.stderr)
    if "lucaprot_also_identified" in out.columns:
        n_also = (out["lucaprot_also_identified"] == "True").sum()
        print(f"  contigs also identified by LucaProt (kept or overlap-dropped): {n_also}",
              file=sys.stderr)


if __name__ == "__main__":
    main()
