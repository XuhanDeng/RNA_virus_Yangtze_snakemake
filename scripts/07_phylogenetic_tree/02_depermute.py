#!/usr/bin/env python3
"""
De-permute RdRp sequences from palm-confirmed FASTA files.

For sequences with permuted C-A-B motif order (C_pos < A_pos), cuts at A_pos
and moves the N-terminal C-containing block to the end, producing canonical A-B-C order.

Input FASTAs must have palm_annot headers:
  >{sequence_id} A:{pos}:{motif} B:{pos}:{motif} C:{pos}:{motif}
  where positions are 1-based.

All three categories (High-confident, RDRPCatch_palm, lucaprot_palmscan) use this
same header format — no TSV lookup needed.

Usage:
  python 02_depermute.py \
    --contig-dir  result/07_phylogenetic_tree/0_input/Contig \
    --ictv-dir    result/07_phylogenetic_tree/0_input/ICTV \
    --out-dir     result/07_phylogenetic_tree/1_depermuted \
    --categories  High-confident RDRPCatch_palm lucaprot_palmscan
"""

import argparse
import os
import re
import sys

LINE_WIDTH = 60

CATEGORIES = ["High-confident", "RDRPCatch_palm", "lucaprot_palmscan"]

# Regex to parse A/B/C from palm_annot header (all three required together)
_ABC_RE = re.compile(
    r"A:(\d+):(\S+)\s+B:(\d+):(\S+)\s+C:(\d+):(\S+)"
)

# Individual motif regexes for partial-hit summary rows
_A_RE = re.compile(r"A:(\d+):(\S+)")
_B_RE = re.compile(r"B:(\d+):(\S+)")
_C_RE = re.compile(r"C:(\d+):(\S+)")


def parse_header(header: str):
    """
    Parse palm_annot header line (without leading '>').
    Returns (seq_id, A_pos_1based, A_motif, B_pos_1based, B_motif, C_pos_1based, C_motif)
    or raises ValueError if not all three motifs are present.
    """
    parts = header.split(None, 1)
    seq_id = parts[0]
    annotation = parts[1] if len(parts) > 1 else ""
    m = _ABC_RE.search(annotation)
    if not m:
        raise ValueError(f"Cannot parse A/B/C from header: {header!r}")
    a_pos, a_mot, b_pos, b_mot, c_pos, c_mot = m.groups()
    return seq_id, int(a_pos), a_mot, int(b_pos), b_mot, int(c_pos), c_mot


def parse_partial_header(header: str):
    """Extract whatever A/B/C motifs are present; return NA for missing ones."""
    parts = header.split(None, 1)
    seq_id = parts[0]
    annotation = parts[1] if len(parts) > 1 else ""
    ma = _A_RE.search(annotation)
    mb = _B_RE.search(annotation)
    mc = _C_RE.search(annotation)
    return (
        seq_id,
        int(ma.group(1)) if ma else "NA", ma.group(2) if ma else "NA",
        int(mb.group(1)) if mb else "NA", mb.group(2) if mb else "NA",
        int(mc.group(1)) if mc else "NA", mc.group(2) if mc else "NA",
    )


def wrap(seq: str, width: int = LINE_WIDTH) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def process_fasta(fasta_path: str, source: str, category: str,
                  out_fh, summary_rows: list):
    """
    Read one input FASTA, de-permute where needed, write to out_fh,
    append summary rows.
    """
    processed = 0
    header = None
    raw_lines = []

    def flush():
        nonlocal processed
        if header is None:
            return
        seq = "".join(line.strip() for line in raw_lines)
        if not seq:
            return

        try:
            seq_id, a1, a_mot, b1, b_mot, c1, c_mot = parse_header(header)
            # Convert 1-based to 0-based
            a_pos = a1 - 1
            if c1 < a1:
                # Permuted C-A-B: move C block to end
                seq_out = seq[a_pos:] + seq[:a_pos]
                depermuted = "yes"
            else:
                seq_out = seq
                depermuted = "no"
            # Sanity check
            for motif, name in [(a_mot, "A"), (b_mot, "B"), (c_mot, "C")]:
                if motif not in seq_out:
                    print(f"  WARNING: {seq_id} motif {name} ({motif}) not found "
                          f"after de-permutation", file=sys.stderr)
        except ValueError:
            # Incomplete motifs (e.g. only B+C or only A+B) — write as-is, no de-permutation
            seq_id, a1, a_mot, b1, b_mot, c1, c_mot = parse_partial_header(header)
            seq_out = seq
            depermuted = "no"

        out_fh.write(f">{seq_id}\n{wrap(seq_out)}\n")

        summary_rows.append({
            "seq_id":      seq_id,
            "source":      source,
            "category":    category,
            "A_pos":       a1,
            "A_motif":     a_mot,
            "B_pos":       b1,
            "B_motif":     b_mot,
            "C_pos":       c1,
            "C_motif":     c_mot,
            "depermuted":  depermuted,
        })
        processed += 1

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                header = line[1:]
                raw_lines = []
            else:
                raw_lines.append(line)
    flush()

    print(f"  {source}/{category}: {processed} sequences", file=sys.stderr)
    return processed


def write_summary(summary_rows: list, path: str):
    cols = ["seq_id", "source", "category",
            "A_pos", "A_motif", "B_pos", "B_motif", "C_pos", "C_motif",
            "depermuted"]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in summary_rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    print(f"Summary: {len(summary_rows)} sequences -> {path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--contig-dir", required=True,
                        help="Input dir with Contig FASTAs (0_input/Contig/)")
    parser.add_argument("--ictv-dir",   required=True,
                        help="Input dir with ICTV FASTAs (0_input/ICTV/)")
    parser.add_argument("--out-dir",    required=True,
                        help="Output root (1_depermuted/)")
    parser.add_argument("--categories", nargs="+", default=CATEGORIES,
                        help="Categories to process")
    args = parser.parse_args()

    contig_out = os.path.join(args.out_dir, "Contig")
    ictv_out   = os.path.join(args.out_dir, "ICTV")
    os.makedirs(contig_out, exist_ok=True)
    os.makedirs(ictv_out,   exist_ok=True)

    summary_rows = []
    merged_path  = os.path.join(args.out_dir, "Contig_ICTV_rdrp.faa")

    with open(merged_path, "w") as merged_fh:
        for source, in_dir, out_subdir in [
            ("Contig", args.contig_dir, contig_out),
            ("ICTV",   args.ictv_dir,   ictv_out),
        ]:
            for cat in args.categories:
                fasta_in  = os.path.join(in_dir,      f"{cat}.faa")
                fasta_out = os.path.join(out_subdir,  f"{cat}.faa")

                if not os.path.exists(fasta_in):
                    print(f"  SKIP (not found): {fasta_in}", file=sys.stderr)
                    continue

                with open(fasta_out, "w") as cat_fh:
                    # Write to both per-category file and merged file
                    class Tee:
                        def write(self, s):
                            cat_fh.write(s)
                            merged_fh.write(s)

                    process_fasta(fasta_in, source, cat, Tee(), summary_rows)

    summary_path = os.path.join(args.out_dir, "Contig_ICTV_rdrp_summary.tsv")
    write_summary(summary_rows, summary_path)
    print(f"Merged FASTA: {merged_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
