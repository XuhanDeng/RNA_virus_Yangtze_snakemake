#!/usr/bin/env python3
"""
Merge Yangtze + ICTV de-permuted sequences with RVMT reference RdRp sequences.

RVMT headers are prefixed with 'RVMT|' to avoid ID collision with Yangtze/ICTV IDs.
Output is Yangtze/ICTV sequences first, then RVMT references.

Usage:
  python 03_add_reference.py \
    --yangtze result/07_phylogenetic_tree/1_depermuted/Contig_ICTV_rdrp.faa \
    --rvmt    database/RVMT/RVMT_Zenodo_V4/RVMT_RdRp.faa \
    --out     result/07_phylogenetic_tree/2_with_ref/Contig_ICTV_RVMT_rdrp.faa
"""

import argparse
import os
import sys

LINE_WIDTH = 60


def copy_fasta(src_path: str, out_fh, prefix: str = ""):
    """Copy a FASTA file to out_fh, optionally prefixing each sequence ID."""
    n = 0
    header = None
    seq_lines = []

    def flush():
        nonlocal n
        if header is None:
            return
        seq = "".join(line.strip() for line in seq_lines)
        if not seq:
            return
        seq_id = header.split()[0]
        rest   = header[len(seq_id):]          # trailing annotation if any
        new_id = prefix + seq_id if prefix else seq_id
        out_fh.write(f">{new_id}{rest}\n")
        for i in range(0, len(seq), LINE_WIDTH):
            out_fh.write(seq[i:i+LINE_WIDTH] + "\n")
        n += 1

    with open(src_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header    = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
    flush()
    return n


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--yangtze", required=True,
                        help="Merged Yangtze+ICTV de-permuted FASTA")
    parser.add_argument("--rvmt",    required=True,
                        help="RVMT RdRp reference FASTA")
    parser.add_argument("--out",     required=True,
                        help="Output merged FASTA path")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    with open(args.out, "w") as out_fh:
        n_yangtze = copy_fasta(args.yangtze, out_fh, prefix="")
        print(f"Yangtze+ICTV sequences: {n_yangtze}", file=sys.stderr)

        n_rvmt = copy_fasta(args.rvmt, out_fh, prefix="RVMT|")
        print(f"RVMT reference sequences: {n_rvmt}", file=sys.stderr)

    print(f"Total: {n_yangtze + n_rvmt} -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
