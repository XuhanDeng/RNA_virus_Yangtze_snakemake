#!/usr/bin/env python3
"""
Extract palm core and extended regions from RdRpCATCH / LucaProt sequences.

Motif info (A/B/C sequences) is read from the palm_annot headers in:
  8_RdRp_protein/RdRpCATCH.faa
  8_RdRp_protein/LucaProt.faa

Full sequences to search motifs in are read from:
  RdRpCATCH: 1_RdRpCATCH/{sample}/*_trimmed_aminoacid_sequences.fasta  (all samples)
  LucaProt:  3_lucaprot/{sample}/*_lucaprot_proteins.faa               (all samples)

Strategy:
  - Parse motif sequences (A/B/C) from the palm_annot header
  - String-search each motif in the full sequence (ignore header positions)
  - Take the leftmost occurrence of each motif
  - Palm core   = min(motif starts) .. max(motif ends)
  - Extended    = core +/- flank aa (capped at sequence boundaries)

Outputs per source (RdRpCATCH / LucaProt):
  {outdir}/{source}/rdrp_full.faa        full sequence, header unchanged
  {outdir}/{source}/palm_core.faa        tight palm core only
  {outdir}/{source}/palm_extended.faa    palm core +/- flank aa
  {outdir}/{source}/palm_regions.tsv     per-sequence table
"""

import argparse
import glob
import re
import sys
from pathlib import Path

FLANK      = 150
MOTIF_KEYS = ["A", "B", "C"]


# ── fasta helpers ──────────────────────────────────────────────────────────────

def read_fasta(path: str) -> dict:
    """Return {seq_id: sequence} from a fasta file. seq_id is everything after '>'
    up to the first space, with 'lcl|' prefix stripped."""
    seqs = {}
    sid = seq = None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if sid is not None:
                    seqs[sid] = seq
                raw = line[1:].split()[0]
                sid = raw.removeprefix("lcl|")
                seq = ""
            elif sid is not None:
                seq += line
    if sid is not None:
        seqs[sid] = seq
    return seqs


def load_fasta_glob(pattern: str) -> dict:
    """Merge all fastas matching glob pattern into one {seq_id: seq} dict."""
    merged = {}
    for path in glob.glob(pattern, recursive=True):
        merged.update(read_fasta(path))
    return merged


def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


# ── motif helpers ──────────────────────────────────────────────────────────────

def parse_motifs_from_header(header: str) -> dict:
    """Return {'A': seq, 'B': seq, 'C': seq} from palm_annot header."""
    motifs = {}
    for key in MOTIF_KEYS:
        m = re.search(rf'\b{key}:(\d+):(\S+)', header)
        motifs[key] = m.group(2) if m else ""
    return motifs


def find_leftmost(seq: str, pattern: str):
    """Return (start, end) 0-based of leftmost hit, or None."""
    if not pattern:
        return None
    idx = seq.find(pattern)
    return (idx, idx + len(pattern)) if idx != -1 else None


def count_occurrences(seq: str, pattern: str) -> int:
    if not pattern:
        return 0
    n, start = 0, 0
    while True:
        idx = seq.find(pattern, start)
        if idx == -1:
            break
        n += 1
        start = idx + 1
    return n


# ── TSV row builder ────────────────────────────────────────────────────────────

TSV_COLS = [
    "seq_id", "source", "rdrp_length",
    "motif_a_seq", "motif_a_start", "motif_a_end", "motif_a_occurrences",
    "motif_b_seq", "motif_b_start", "motif_b_end", "motif_b_occurrences",
    "motif_c_seq", "motif_c_start", "motif_c_end", "motif_c_occurrences",
    "core_start", "core_end", "core_length",
    "extended_start", "extended_end", "extended_length",
    "flag",
]


def make_row(seq_id, source, full_seq, motifs, hits, occ,
             core_start, core_end, ext_start, ext_end, flags):
    def pos(hit):
        return (str(hit[0]), str(hit[1])) if hit else ("", "")

    a0, a1 = pos(hits.get("A"))
    b0, b1 = pos(hits.get("B"))
    c0, c1 = pos(hits.get("C"))

    core_len = (core_end - core_start) if core_start is not None else ""
    ext_len  = (ext_end  - ext_start)  if ext_start  is not None else ""

    return "\t".join([
        seq_id, source, str(len(full_seq)),
        motifs.get("A", ""), a0, a1, str(occ.get("A", 0)),
        motifs.get("B", ""), b0, b1, str(occ.get("B", 0)),
        motifs.get("C", ""), c0, c1, str(occ.get("C", 0)),
        str(core_start) if core_start is not None else "",
        str(core_end)   if core_end   is not None else "",
        str(core_len),
        str(ext_start)  if ext_start  is not None else "",
        str(ext_end)    if ext_end    is not None else "",
        str(ext_len),
        ";".join(flags) if flags else "ok",
    ])


# ── main processing ────────────────────────────────────────────────────────────

def process(motif_fasta: str, full_seq_db: dict, outdir: Path,
            source: str, flank: int):
    outdir.mkdir(parents=True, exist_ok=True)

    f_full = open(outdir / "rdrp_full.faa",     "w")
    f_core = open(outdir / "palm_core.faa",     "w")
    f_ext  = open(outdir / "palm_extended.faa", "w")
    f_tsv  = open(outdir / "palm_regions.tsv",  "w")
    f_tsv.write("\t".join(TSV_COLS) + "\n")

    total = ok = flagged = missing = 0

    with open(motif_fasta) as fh:
        header = full_header = seq_lines = None

        def emit(full_header, motif_seq_unused):
            nonlocal ok, flagged, missing, total
            total += 1

            seq_id = full_header.split()[0]
            motifs = parse_motifs_from_header(full_header)

            # look up full sequence
            full_seq = full_seq_db.get(seq_id)
            if full_seq is None:
                missing += 1
                print(f"  MISSING full seq: {seq_id}", file=sys.stderr)
                return

            hits  = {}
            occ   = {}
            flags = []

            for key in MOTIF_KEYS:
                mseq     = motifs[key]
                hit      = find_leftmost(full_seq, mseq)
                n        = count_occurrences(full_seq, mseq)
                hits[key] = hit
                occ[key]  = n
                if not mseq:
                    flags.append(f"{key}_missing_in_header")
                elif hit is None:
                    flags.append(f"{key}_not_found_in_seq")
                elif n > 1:
                    flags.append(f"{key}_duplicate({n})")

            found = {k: v for k, v in hits.items() if v is not None}

            # always write full sequence
            f_full.write(f">{full_header}\n{wrap(full_seq)}\n")

            if not found:
                flags.append("no_motifs_found_skip_trim")
                flagged += 1
                f_core.write(f">{full_header}\n{wrap(full_seq)}\n")
                f_ext.write(f">{full_header}\n{wrap(full_seq)}\n")
                f_tsv.write(make_row(seq_id, source, full_seq, motifs,
                                     hits, occ, None, None, None, None, flags) + "\n")
                return

            core_start = min(v[0] for v in found.values())
            core_end   = max(v[1] for v in found.values())
            ext_start  = max(0, core_start - flank)
            ext_end    = min(len(full_seq), core_end + flank)

            f_core.write(f">{full_header}\n{wrap(full_seq[core_start:core_end])}\n")
            f_ext.write(f">{full_header}\n{wrap(full_seq[ext_start:ext_end])}\n")

            if flags:
                flagged += 1
            else:
                ok += 1

            f_tsv.write(make_row(seq_id, source, full_seq, motifs,
                                  hits, occ, core_start, core_end,
                                  ext_start, ext_end, flags) + "\n")

        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if full_header is not None:
                    emit(full_header, "".join(seq_lines))
                full_header = line[1:]
                seq_lines   = []
            elif full_header is not None:
                seq_lines.append(line)
        if full_header is not None:
            emit(full_header, "".join(seq_lines))

    for f in (f_full, f_core, f_ext, f_tsv):
        f.close()

    print(f"[{source}] total={total}  ok={ok}  flagged={flagged}  "
          f"missing_full_seq={missing}", file=sys.stderr)


# ── CLI ────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rdrpcatch-motif",
                        help="8_RdRp_protein/RdRpCATCH.faa (motif headers)")
    parser.add_argument("--rdrpcatch-seqs",
                        help="Glob to 1_RdRpCATCH/**/*_trimmed_aminoacid_sequences.fasta")
    parser.add_argument("--lucaprot-motif",
                        help="8_RdRp_protein/LucaProt.faa (motif headers)")
    parser.add_argument("--lucaprot-seqs",
                        help="Glob to 3_lucaprot/**/*_lucaprot_proteins.faa")
    parser.add_argument("--outdir",  required=True)
    parser.add_argument("--flank",   type=int, default=150)
    args = parser.parse_args()

    outdir = Path(args.outdir)

    if args.rdrpcatch_motif:
        if not args.rdrpcatch_seqs:
            parser.error("--rdrpcatch-seqs required when --rdrpcatch-motif given")
        print("Loading RdRpCATCH full sequences...", file=sys.stderr)
        db = load_fasta_glob(args.rdrpcatch_seqs)
        print(f"  {len(db)} sequences loaded", file=sys.stderr)
        process(args.rdrpcatch_motif, db, outdir / "RdRpCATCH", "RdRpCATCH", args.flank)

    if args.lucaprot_motif:
        if not args.lucaprot_seqs:
            parser.error("--lucaprot-seqs required when --lucaprot-motif given")
        print("Loading LucaProt full sequences...", file=sys.stderr)
        db = load_fasta_glob(args.lucaprot_seqs)
        print(f"  {len(db)} sequences loaded", file=sys.stderr)
        process(args.lucaprot_motif, db, outdir / "LucaProt", "LucaProt", args.flank)


if __name__ == "__main__":
    main()
