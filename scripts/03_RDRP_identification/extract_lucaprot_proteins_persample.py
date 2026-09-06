#!/usr/bin/env python3
"""
Write a FASTA from a single per-sample lucaprot CSV for palm_annot input.

The lucaprot CSV (RdRPs_only_using_threshold*.csv) has columns:
  protein_id, seq, prob, label

protein_id example:
  >lcl|ORF9_AnQing_D_33_0000000002:8171:2670 unnamed protein product

FASTA header written: >ORF9_AnQing_D_33_0000000002:8171:2670
(stripped of ">lcl|" prefix and description text)

This FASTA is used as input to palm_annot (palm_annot_lucaprot_persample rule).
The resulting palmscan TSV Label column will match these headers, allowing
merge_rdrp_results.py to link palmscan hits back to the original ORFs.
"""

import argparse
import re
import sys

import pandas as pd


def clean_header(protein_id: str) -> str:
    s = str(protein_id).strip().lstrip(">")
    s = re.sub(r"^lcl\|", "", s)
    return s.split()[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv",    required=True, help="RdRPs_only_using_threshold*.csv")
    parser.add_argument("--output", required=True, help="Output FASTA")
    args = parser.parse_args()

    df = pd.read_csv(args.csv, dtype=str)
    if df.empty or "protein_id" not in df.columns or "seq" not in df.columns:
        open(args.output, "w").close()
        print("Empty or missing columns — wrote empty FASTA", file=sys.stderr)
        return

    written = 0
    with open(args.output, "w") as fh:
        for _, row in df.iterrows():
            header = clean_header(row["protein_id"])
            seq = str(row["seq"]).strip()
            if seq and seq != "nan":
                fh.write(f">{header}\n{seq}\n")
                written += 1

    print(f"Written: {written} sequences -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
