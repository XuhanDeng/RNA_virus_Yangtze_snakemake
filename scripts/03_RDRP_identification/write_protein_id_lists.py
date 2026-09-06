#!/usr/bin/env python3
"""
Write per-source protein ID list files from the combined contig table.

Two ID files, used by seqkit grep to pull sequences out of the matching
source FASTA:
  RdRpCATCH_ids.txt — ID format {contig}_frame={N}_RdRp_{start}-{end}
                       extracted from rdrp_trimmed.faa (palmscan output on RC)
  LucaProt_ids.txt   — ID format ORF{n}_{contig}:{start}:{end}
                       extracted from lucaprot_rdrp_trimmed.faa (palmscan output on LP)
"""

import argparse
import os
import re
import sys

import pandas as pd

ALL_ID_FILES = ["RdRpCATCH", "LucaProt"]


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
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged", required=True,
                        help="all_samples_rdrp_merged.tsv (or ICTV_rdrp_merged.tsv)")
    parser.add_argument("--outdir", required=True,
                        help="Directory to write {source}_ids.txt files")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_csv(args.merged, sep="\t", dtype=str)
    print(f"Loaded {len(df)} contigs from {args.merged}", file=sys.stderr)

    ids = {"RdRpCATCH": set(), "LucaProt": set()}
    for _, row in df.iterrows():
        source = str(row.get("source", "")).strip()
        if source == "RdRpCATCH":
            pid = _rc_id(row)
            if pid:
                ids["RdRpCATCH"].add(pid)
        elif source == "LucaProt":
            pid = _lp_id(row)
            if pid:
                ids["LucaProt"].add(pid)

    for key in ALL_ID_FILES:
        path = os.path.join(args.outdir, f"{key}_ids.txt")
        sorted_ids = sorted(ids[key])
        with open(path, "w") as fh:
            fh.write("\n".join(sorted_ids) + ("\n" if sorted_ids else ""))
        print(f"  {key}: {len(sorted_ids)} IDs -> {path}", file=sys.stderr)


if __name__ == "__main__":
    main()
