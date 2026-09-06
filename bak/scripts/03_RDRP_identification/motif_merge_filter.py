#!/usr/bin/env python3
"""
motif_merge_filter.py
Merge, filter and classify RdRP motif hits from hmmsearch, PSI-BLAST, MMseqs2, Diamond.

Usage:
    python motif_merge_filter.py --search_dir search/ --input_faa all_candidates.faa --outdir motif_results/

Input structure:
    search/
    ├── hmmsearch/  mot.1.tsv  mot.2.tsv  mot.3.tsv  mot.4.tsv
    ├── psiblast/   mot.1.tsv  mot.2.tsv  mot.3.tsv  mot.4.tsv
    ├── mmseqs/     mot.1.tsv  mot.2.tsv  mot.3.tsv  mot.4.tsv
    └── diamond/    mot.1.tsv  mot.2.tsv  mot.3.tsv  mot.4.tsv
"""

import argparse
import os
import re
import sys
import pandas as pd
from Bio import SeqIO

# ─────────────────────────────────────────────────────────────────────────────
# Search tool orientation for this pipeline:
#   hmmsearch : motif HMM      → query ;  your seq → target
#   psiblast  : motif MSA      → query ;  your seq → subject
#   mmseqs    : motif profile  → query ;  your seq → target
#   diamond   : motif seq      → query ;  your seq → target
#
# All parsers extract coordinates so that in the output:
#   q1/q2 = position on YOUR sequence
#   p1/p2 = position on motif profile/sequence
# ─────────────────────────────────────────────────────────────────────────────


def parse_hmmsearch(path):
    # Input TSV pattern (domtblout, space-separated, comment lines start with #):
    #   [0]          [1]   [2]  [3]         [4]   [5]  [6]     [7]    [8]   [9]  [10] [11]      [12]     [13]     [14]       [15] [16] [17] [18]     [19]     [20]   [21]  [22]
    #   target_name  tacc  qL   query_name  qacc  pL   evalue  score  bias  #    of   c-evalue  i-evalue score2   bias2      p1   p2   q1   q2       env_from env_to acc   description
    #   (your seq)              (motif)
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.split()
            if len(parts) < 19:
                continue
            rows.append({
                "RdRp_id": parts[0],    # target_name = YOUR sequence
                "profile": parts[3],    # query_name  = motif
                "qL":      parts[2],
                "pL":      parts[5],
                "evalue":  parts[6],
                "score":   parts[7],
                "p1":      parts[15],   # hmm coord from
                "p2":      parts[16],   # hmm coord to
                "q1":      parts[17],   # ali coord on YOUR seq
                "q2":      parts[18],   # ali coord on YOUR seq
            })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def parse_psiblast(path):
    # Input TSV pattern (tab-separated, outfmt6 with sseqid first + profile appended as col 12):
    #   [0]      [1]     [2]   [3]   [4]   [5]   [6]   [7]   [8]   [9]      [10]   [11]
    #   sseqid   pident  q1    q2    p1    p2    qL    pL    ali   evalue   score  profile
    #   (your seq)       (q1)  (q2)  (p1)  (p2)  (qL)  (pL)  (ali)         (bits)  (motif name)
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 12:
                continue
            rows.append({
                "RdRp_id": parts[0],   # sseqid = YOUR sequence
                "q1":      parts[2],   # sstart = position on YOUR seq
                "q2":      parts[3],   # send   = position on YOUR seq
                "p1":      parts[4],   # qstart = position on motif
                "p2":      parts[5],   # qend   = position on motif
                "qL":      parts[6],   # slen
                "pL":      parts[7],   # qlen
                "ali_len": parts[8],
                "evalue":  parts[9],
                "score":   parts[10],
                "profile": parts[11],
            })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def parse_mmseqs(path):
    # Input TSV pattern (tab-separated, mmseqs convertalis format):
    #   [0]    [1]     [2]     [3]      [4]     [5]     [6]     [7]   [8]   [9]     [10]  [11]  [12]    [13]  [14]  [15]    [16]      [17]  [18]
    #   query  target  evalue  gapopen  pident  nident  qstart  qend  qlen  tstart  tend  tlen  alnlen  raw   bits  qframe  mismatch  qcov  tcov
    #   (motif)  (your seq)              (p1)    (p2)  (pL)  (q1)   (q2)  (qL)  (ali)         (score)
    #
    # YOUR setup: query = motif profile, target = YOUR RdRP sequence
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            rows.append({
                "profile": parts[0],    # query  = motif
                "RdRp_id": parts[1],    # target = YOUR sequence
                "evalue":  parts[2],
                "p1":      parts[6],    # qstart = position on motif
                "p2":      parts[7],    # qend   = position on motif
                "pL":      parts[8],    # qlen
                "q1":      parts[9],    # tstart = position on YOUR seq
                "q2":      parts[10],   # tend   = position on YOUR seq
                "qL":      parts[11],   # tlen
                "ali_len": parts[12],
                "score":   parts[14],   # bits
            })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


def parse_diamond(path):
    # Input TSV pattern (tab-separated, diamond blastp outfmt6):
    #   [0]     [1]     [2]     [3]      [4]     [5]     [6]     [7]   [8]     [9]   [10]  [11]  [12]      [13]
    #   qseqid  sseqid  evalue  gapopen  pident  length  qstart  qend  sstart  send  qlen  slen  mismatch  bitscore
    #   (motif) (your seq)               (ali)   (p1)    (p2)   (q1)  (q2)   (pL)  (qL)            (score)
    #
    # YOUR setup: query = motif sequence, target = YOUR RdRP sequence
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or line.strip() == "":
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 14:
                continue
            rows.append({
                "profile": parts[0],    # qseqid = motif
                "RdRp_id": parts[1],    # sseqid = YOUR sequence
                "evalue":  parts[2],
                "ali_len": parts[5],    # length = alignment length
                "p1":      parts[6],    # qstart = position on motif
                "p2":      parts[7],    # qend   = position on motif
                "q1":      parts[8],    # sstart = position on YOUR seq
                "q2":      parts[9],    # send   = position on YOUR seq
                "pL":      parts[10],   # qlen   = motif length
                "qL":      parts[11],   # slen   = your seq length
                "score":   parts[13],   # bitscore
            })
    return pd.DataFrame(rows) if rows else pd.DataFrame()


PARSERS = {
    "hmmsearch": parse_hmmsearch,
    "psiblast":  parse_psiblast,
    "mmseqs":    parse_mmseqs,
    "diamond":   parse_diamond,
}

# ─────────────────────────────────────────────────────────────────────────────
# All parsers return a DataFrame with these standardised columns:
#
#   RdRp_id | profile | q1 | q2 | p1 | p2 | qL | pL | score | evalue | ali_len
#   ───────   ───────   ──   ──   ──   ──   ──   ──   ─────   ──────   ───────
#   your seq  motif     start/end on  start/end   your  motif  bit      E-      alignment
#   ID        cluster   YOUR sequence on motif    seq   len    score    value   length
#                       (AA position) profile     len
#
# All 3 tools run motif as query, YOUR sequence as target/subject:
#   hmmsearch : motif HMM  → query ;  your seq → target  (domtblout q1/q2 = ali coords on target)
#   psiblast  : motif MSA  → query ;  your seq → subject (sstart/send = coords on your seq)
#   mmseqs    : motif prof → query ;  your seq → target  (tstart/tend = coords on your seq)
#
# q1/q2 in output always = position on YOUR sequence
# p1/p2 in output always = position on motif profile
# ─────────────────────────────────────────────────────────────────────────────


def load_all_hits(search_dir, tools, motifs):
    frames = []
    for tool in tools:
        for m in motifs:
            path = os.path.join(search_dir, tool, f"mot.{m}.tsv")
            if not os.path.exists(path):
                print(f"  [skip] missing: {path}")
                continue
            if os.path.getsize(path) == 0:
                print(f"  [skip] empty:   {path}")
                continue
            df = PARSERS[tool](path)
            if df is None or df.empty:
                print(f"  [skip] no hits: {path}")
                continue
            df["motif_type"] = m
            df["search_tool"] = tool
            frames.append(df)
            print(f"  loaded {len(df):>6} hits  {tool}/mot.{m}.tsv")
    if not frames:
        sys.exit("ERROR: no hits loaded from any file")
    return pd.concat(frames, ignore_index=True)


def coerce_types(df):
    for col in ["evalue", "score", "pident"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    for col in ["q1", "q2", "p1", "p2", "qL", "pL", "ali_len"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").astype("Int64")
    df["ali_len"] = (df["q2"] - df["q1"]).abs()
    return df


def filter_hits(df):
    # Neri thresholds — Annotations.r lines 907-909
    n0 = len(df)
    df = df[df["evalue"]  < 1e-5]
    df = df[df["score"]   > 10  ]
    df = df[df["ali_len"] > 10  ]
    print(f"  {len(df)} hits kept  (from {n0} raw)")
    return df


def best_hit_per_motif(df):
    df = df.sort_values(["RdRp_id", "motif_type", "score", "evalue"],
                        ascending=[True, True, False, True])
    best = df.groupby(["RdRp_id", "motif_type"], as_index=False).first()
    print(f"  {len(best)} best hits  ({best['RdRp_id'].nunique()} unique sequences)")
    return best


def get_contig_name(rid):
    """Extract base contig name from either rdrpcatch or lucaprot ID."""
    if rid.startswith("ORF"):
        # ORFn_ContigName:start:end  →  ContigName
        return re.sub(r"^ORF\d+_", "", rid).split(":")[0]
    else:
        # ContigName_frame=X  →  ContigName
        return re.sub(r"_frame=.*", "", rid)


def is_lucaprot(rid):
    return rid.startswith("ORF")


def deduplicate_by_contig(hits_best):
    """
    For each contig, if both rdrpcatch and lucaprot sequences exist,
    keep only the rdrpcatch sequence. Motif coordinates are protein-internal
    and not comparable across different ORF calls, so same-contig presence
    of both is sufficient evidence they represent the same RdRp.
    """
    n_before = hits_best["RdRp_id"].nunique()

    hits_best = hits_best.copy()
    hits_best["_contig"] = hits_best["RdRp_id"].apply(get_contig_name)
    hits_best["_is_lp"]  = hits_best["RdRp_id"].apply(is_lucaprot)

    # Contigs that have at least one rdrpcatch hit
    rc_contigs = set(hits_best.loc[~hits_best["_is_lp"], "_contig"])

    # Drop all lucaprot sequences whose contig is already covered by rdrpcatch
    drop_ids = set(
        hits_best.loc[
            hits_best["_is_lp"] & hits_best["_contig"].isin(rc_contigs),
            "RdRp_id"
        ]
    )

    hits_dedup = hits_best[~hits_best["RdRp_id"].isin(drop_ids)].drop(
        columns=["_contig", "_is_lp"]
    )
    n_after = hits_dedup["RdRp_id"].nunique()
    print(f"  Dropped {len(drop_ids)} lucaprot sequences from contigs covered by rdrpcatch")
    print(f"  {n_before} → {n_after} unique sequences after deduplication")
    return hits_dedup


def extract_motif_seqs(hits_best, input_faa):
    seqs = {r.id: r.seq for r in SeqIO.parse(input_faa, "fasta")}
    records = []
    for _, row in hits_best.iterrows():
        seq_id  = row["RdRp_id"]
        profile = row["profile"]
        q1      = int(row["q1"]) - 1   # convert to 0-based
        q2      = int(row["q2"])
        if seq_id not in seqs:
            continue
        subseq = seqs[seq_id][q1:q2]
        name   = f"{seq_id}.{profile}"
        records.append((name, str(subseq)))
    return records


def classify_motif_order(hits_best):
    # Pivot: one row per RdRp_id, q1 position for each motif
    pivot = hits_best.pivot_table(
        index="RdRp_id",
        columns="motif_type",
        values="q1",
        aggfunc="min"
    )
    pivot.columns = [f"q1_{['A','B','C','D'][int(c)-1]}" for c in pivot.columns]
    pivot = pivot.reset_index()

    # Attach contig and source (take first value per RdRp_id — same for all rows)
    meta = hits_best.groupby("RdRp_id")[["contig", "source"]].first().reset_index()
    pivot = pivot.merge(meta, on="RdRp_id", how="left")

    # Count detected motifs
    pos_cols = [c for c in ["q1_A","q1_B","q1_C","q1_D"] if c in pivot.columns]
    pivot["n_motifs"] = pivot[pos_cols].notna().sum(axis=1)

    # Canonical: all 4 present AND A < B < C < D
    if all(c in pivot.columns for c in ["q1_A","q1_B","q1_C","q1_D"]):
        pivot["canonical"] = (
            pivot["q1_A"].notna() & pivot["q1_B"].notna() &
            pivot["q1_C"].notna() & pivot["q1_D"].notna() &
            (pivot["q1_A"] < pivot["q1_B"]) &
            (pivot["q1_B"] < pivot["q1_C"]) &
            (pivot["q1_C"] < pivot["q1_D"])
        )
        # Permuted: motif C upstream of motif A  (C-A-B-D)
        pivot["permuted"] = (
            pivot["q1_C"].notna() & pivot["q1_A"].notna() &
            (pivot["q1_C"] < pivot["q1_A"])
        )
    else:
        pivot["canonical"] = None
        pivot["permuted"]  = None

    pivot["incomplete"] = pivot["n_motifs"] < 4

    # motif_order: letters of detected motifs sorted by q1 position, e.g. "ABCD", "CABD", "AB"
    letter_map = {"q1_A": "A", "q1_B": "B", "q1_C": "C", "q1_D": "D"}
    def _order_string(row):
        present = {col: row[col] for col in letter_map if col in row.index and pd.notna(row[col])}
        return "".join(letter_map[col] for col in sorted(present, key=lambda c: present[c]))
    pivot["motif_order"] = pivot.apply(_order_string, axis=1)

    return pivot


def depermute_sequences(hits_best, motif_order, input_faa):
    # For permuted RdRPs (C-A-B-D order), rearrange the full-length sequence so
    # motifs appear in A-B-C-D order. The permuted arrangement has the C-D block
    # upstream of the A-B block, so we swap those two halves:
    #
    #   Original:  [C...D | A...B]   (positions on protein)
    #   Output:    [A...B | C...D]   (depermuted, now A-B-C-D order)
    #
    # The cut point is the midpoint between the end of D and the start of A,
    # i.e. between the two blocks. We use q2 of D and q1 of A (0-based).
    #
    # Canonical sequences are written unchanged.
    # Input motif_order must already be filtered to complete sequences only.

    seqs = {r.id: str(r.seq) for r in SeqIO.parse(input_faa, "fasta")}

    # Build lookup: RdRp_id -> q1/q2 per motif
    motif_coords = {}
    for _, row in hits_best.iterrows():
        rid = row["RdRp_id"]
        m   = int(row["motif_type"])
        if rid not in motif_coords:
            motif_coords[rid] = {}
        motif_coords[rid][m] = (int(row["q1"]), int(row["q2"]))

    records = []
    n_canonical = n_permuted = n_skipped = 0

    for _, row in motif_order.iterrows():
        rid = row["RdRp_id"]
        if rid not in seqs:
            continue
        seq = seqs[rid]

        if row.get("permuted", False):
            coords = motif_coords.get(rid, {})
            if not all(m in coords for m in [1, 2, 3, 4]):
                n_skipped += 1
                continue

            q1_A = coords[1][0]
            q1_B = coords[2][0]
            q1_C = coords[3][0]
            q1_D = coords[4][0]

            # Only depermute if D is the last motif
            if not (q1_D > q1_A and q1_D > q1_B and q1_D > q1_C):
                print(f"  [skip] {rid}: D is not last, invalid order")
                n_skipped += 1
                continue

            order = row.get("motif_order", "")

            if order == "CABD":
                # C jumped before A; move C block [q1_C:q1_A-1] to before D
                # Original : [C][A...B...][D...]
                # Depermuted: [A...B...before_D] + [C block] + [D...]
                depermuted = seq[q1_A : q1_D-1] + seq[q1_C : q1_A-1] + seq[q1_D :]

            elif order == "CBAD":
                # C before B, B before A; move both C and B blocks
                # Original : [C][B][A...][D...]
                # Depermuted: [A...before_D] + [B block] + [C block] + [D...]
                depermuted = seq[q1_A : q1_D-1] + seq[q1_B : q1_A-1] + seq[q1_C : q1_B-1] + seq[q1_D :]

            elif order == "ACBD":
                # C jumped between A and B; move C block [q1_C:q1_B-1] to before D
                # Original : [A...][C][B...][D...]
                # Depermuted: [A...before_C] + [B...before_D] + [C block] + [D...]
                depermuted = seq[q1_A : q1_C-1] + seq[q1_B : q1_D-1] + seq[q1_C : q1_B-1] + seq[q1_D :]

            elif order == "BACD":
                # B jumped before A; move B block [q1_B:q1_A-1] to after A
                # Original : [B][A...][C...][D...]
                # Depermuted: [A...before_C] + [B block] + [C...before_D] + [D...]
                depermuted = seq[q1_A : q1_C-1] + seq[q1_B : q1_A-1] + seq[q1_C : q1_D-1] + seq[q1_D :]

            elif order == "BCAD":
                # B and C both before A; move BC block [q1_B:q1_A-1] to before D
                # Original : [B][C][A...][D...]
                # Depermuted: [A...before_D] + [BC block] + [D...]
                depermuted = seq[q1_A : q1_D-1] + seq[q1_B : q1_A-1] + seq[q1_D :]

            else:
                print(f"  [warn] {rid}: unrecognised permutation order '{order}', writing as-is")
                records.append((rid, seq))
                n_canonical += 1
                continue
            records.append((f"{rid}_depermuted", depermuted))
            n_permuted += 1
        else:
            records.append((rid, seq))
            n_canonical += 1

    print(f"  Canonical written as-is : {n_canonical}")
    print(f"  Permuted depermuted     : {n_permuted}")
    print(f"  Incomplete skipped      : {n_skipped}")
    return records


def main():
    parser = argparse.ArgumentParser(description="Merge and filter RdRP motif hits")
    parser.add_argument("--search_dir", required=True,
                        help="Directory with tool subdirs (hmmsearch/, psiblast/, mmseqs/)")
    parser.add_argument("--input_faa",  required=True,
                        help="Full-length candidate RdRP sequences (.faa)")
    parser.add_argument("--outdir",     default="motif_results",
                        help="Output directory")
    parser.add_argument("--tools",      nargs="+",
                        default=["hmmsearch", "psiblast", "mmseqs", "diamond"],
                        help="Tools to include (default: hmmsearch psiblast mmseqs diamond)")
    parser.add_argument("--motifs",     nargs="+", type=int, default=[1,2,3,4],
                        help="Motif numbers to include (default: 1 2 3 4)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ── Step 1: Parse ──
    print("\n=== Step 1: Parsing search results ===")
    all_hits = load_all_hits(args.search_dir, args.tools, args.motifs)
    all_hits = coerce_types(all_hits)
    print(f"  Total raw hits: {len(all_hits)}")

    # ── Step 2: Filter ──
    print("\n=== Step 2: Filtering (evalue<1e-5, score>10, ali_len>10) ===")
    hits_filt = filter_hits(all_hits)

    # ── Step 3: Best hit per sequence per motif ──
    print("\n=== Step 3: Best hit per sequence per motif ===")
    hits_best = best_hit_per_motif(hits_filt)

    # ── Step 3b: Deduplicate rdrpcatch vs lucaprot overlaps ──
    print("\n=== Step 3b: Deduplicating rdrpcatch/lucaprot overlaps ===")
    hits_best = deduplicate_by_contig(hits_best)
    hits_best["contig"] = hits_best["RdRp_id"].apply(get_contig_name)
    hits_best["source"] = hits_best["RdRp_id"].apply(lambda x: "lucaprot" if is_lucaprot(x) else "rdrpcatch")

    out_hits = os.path.join(args.outdir, "motif_hits_best.tsv")
    hits_best.to_csv(out_hits, sep="\t", index=False)
    print(f"  Written: {out_hits}")

    # ── Step 4: Extract motif subsequences ──
    print("\n=== Step 4: Extracting motif subsequences ===")
    records = extract_motif_seqs(hits_best, args.input_faa)
    out_faa = os.path.join(args.outdir, "motif_sequences.faa")
    with open(out_faa, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")
    print(f"  Written: {out_faa}  ({len(records)} sequences)")

    # ── Step 5: Classify canonical / permuted ──
    print("\n=== Step 5: Classifying motif order ===")
    motif_order = classify_motif_order(hits_best)
    out_order = os.path.join(args.outdir, "motif_order.tsv")
    motif_order.to_csv(out_order, sep="\t", index=False)
    print(f"  Written: {out_order}")

    # ── Step 6: Filter to complete sequences (all 4 motifs, D must be last) ──
    print("\n=== Step 6: Filtering complete sequences (all 4 motifs, D last) ===")
    motif_order_complete = motif_order[
        (motif_order["incomplete"] == False) &
        (motif_order["motif_order"].str.endswith("D"))
    ].copy()
    motif_order_suspicious = motif_order[
        (motif_order["incomplete"] == False) &
        (~motif_order["motif_order"].str.endswith("D"))
    ].copy()
    out_order_complete = os.path.join(args.outdir, "motif_order_complete.tsv")
    motif_order_complete.to_csv(out_order_complete, sep="\t", index=False)
    print(f"  {len(motif_order_complete)} complete valid sequences (D last)")
    print(f"  {len(motif_order_suspicious)} suspicious sequences (D not last)")
    print(f"  Written: {out_order_complete}")

    # ── Step 7: Depermute CABD only ──
    print("\n=== Step 7: Depermuting CABD sequences ===")
    depermed = depermute_sequences(hits_best, motif_order_complete, args.input_faa)
    out_depermed = os.path.join(args.outdir, "sequences_depermuted.faa")
    with open(out_depermed, "w") as f:
        for name, seq in depermed:
            f.write(f">{name}\n{seq}\n")
    print(f"  Written: {out_depermed}  ({len(depermed)} sequences)")

    # ── Step 8: Tier 1a — ABCD canonical (original full-length) ──
    print("\n=== Step 8: Tier 1a — ABCD canonical sequences ===")
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(args.input_faa, "fasta")}
    abcd_ids = set(
        motif_order_complete.loc[motif_order_complete["motif_order"] == "ABCD", "RdRp_id"]
    )
    out_tier1a = os.path.join(args.outdir, "tier1a_ABCD.faa")
    n_tier1a = 0
    with open(out_tier1a, "w") as f:
        for rid in abcd_ids:
            if rid in seqs:
                f.write(f">{rid}\n{seqs[rid]}\n")
                n_tier1a += 1
    print(f"  Written: {out_tier1a}  ({n_tier1a} sequences)")

    # ── Step 9: Tier 1b — all depermuted permuted orders (CABD, CBAD, ACBD,
    # BACD, BCAD) — any D-last 4-motif sequence successfully depermuted by
    # Step 7 counts as tier1b, not just the canonical CABD pattern.
    print("\n=== Step 9: Tier 1b — depermuted permuted sequences ===")
    _permuted_orders = {"CABD", "CBAD", "ACBD", "BACD", "BCAD"}
    tier1b_ids = set(
        motif_order_complete.loc[
            motif_order_complete["motif_order"].isin(_permuted_orders), "RdRp_id"
        ]
    )
    out_tier1b = os.path.join(args.outdir, "tier1b_CABD_depermuted.faa")
    n_tier1b = 0
    with open(out_tier1b, "w") as f:
        for name, seq in depermed:
            orig_id = name.replace("_depermuted", "")
            if orig_id in tier1b_ids:
                f.write(f">{name}\n{seq}\n")
                n_tier1b += 1
    print(f"  Written: {out_tier1b}  ({n_tier1b} sequences)")

    # ── Step 10: Tier 2 — 3-motif partial (BCD, ABC, ABD, ACD) ──
    print("\n=== Step 10: Tier 2 — 3-motif partial sequences ===")
    _3motif_orders = {"BCD", "ABC", "ABD", "ACD"}
    tier2_ids = set(
        motif_order.loc[motif_order["motif_order"].isin(_3motif_orders), "RdRp_id"]
    )
    out_tier2 = os.path.join(args.outdir, "tier2_3motif.faa")
    n_tier2 = 0
    with open(out_tier2, "w") as f:
        for rid in tier2_ids:
            if rid in seqs:
                f.write(f">{rid}\n{seqs[rid]}\n")
                n_tier2 += 1
    print(f"  Written: {out_tier2}  ({n_tier2} sequences)")

    # ── Step 11: Suspicious — D not last ──
    print("\n=== Step 11: Suspicious sequences (D not last) ===")
    out_suspicious = os.path.join(args.outdir, "motif_order_suspicious.tsv")
    motif_order_suspicious.to_csv(out_suspicious, sep="\t", index=False)
    print(f"  Written: {out_suspicious}  ({len(motif_order_suspicious)} sequences)")

    # ── Step 12: Tier 3 — all candidates not in tier1/tier2, deduplicated by contig ──
    print("\n=== Step 12: Tier 3 — no-motif candidates ===")
    tier123_ids = abcd_ids | tier1b_ids | tier2_ids

    # Read all sequences from input_faa
    all_input_seqs = list(SeqIO.parse(args.input_faa, "fasta"))

    # Deduplicate all input sequences by contig: keep rdrpcatch, drop lucaprot if same contig covered
    rc_contigs_all = set(
        get_contig_name(r.id) for r in all_input_seqs if not is_lucaprot(r.id)
    )

    out_tier3 = os.path.join(args.outdir, "tier3_no_motif.faa")
    n_tier3 = 0
    with open(out_tier3, "w") as f:
        for r in all_input_seqs:
            rid = r.id
            if rid in tier123_ids:
                continue
            if is_lucaprot(rid) and get_contig_name(rid) in rc_contigs_all:
                continue
            f.write(f">{rid}\n{str(r.seq)}\n")
            n_tier3 += 1
    print(f"  Written: {out_tier3}  ({n_tier3} sequences)")

    # ── Step 13: Comprehensive summary table (tier1+tier2 only) ──
    print("\n=== Step 12: Building comprehensive summary table ===")

    # Assign tier label — tier1/tier2 from motif results, tier3 from input_faa remainder
    tier12_ids = list(abcd_ids) + list(tier1b_ids) + list(tier2_ids - abcd_ids - tier1b_ids)
    tier12_labels = (["tier1a"] * len(abcd_ids) +
                     ["tier1b"] * len(tier1b_ids) +
                     ["tier2"]  * len(tier2_ids - abcd_ids - tier1b_ids))
    tier3_ids = [r.id for r in all_input_seqs
                 if r.id not in tier123_ids
                 and not (is_lucaprot(r.id) and get_contig_name(r.id) in rc_contigs_all)]
    all_ids    = tier12_ids    + tier3_ids
    all_labels = tier12_labels + ["tier3"] * len(tier3_ids)
    tier_ser = pd.DataFrame({"RdRp_id": all_ids, "tier": all_labels})

    # Wide per-motif table from hits_best
    letter = {1: "A", 2: "B", 3: "C", 4: "D"}
    wide_frames = []
    for m, ml in letter.items():
        sub = hits_best[hits_best["motif_type"] == m][
            ["RdRp_id", "q1", "q2", "evalue", "score", "search_tool", "profile"]
        ].copy()
        sub.columns = ["RdRp_id", f"mot{ml}_q1", f"mot{ml}_q2",
                       f"mot{ml}_evalue", f"mot{ml}_score", f"mot{ml}_tool", f"mot{ml}_profile"]
        wide_frames.append(sub.set_index("RdRp_id"))
    wide = pd.concat(wide_frames, axis=1).reset_index()

    # Protein length from input_faa (covers tier3 which has no hits_best entry)
    seq_len_map = {r.id: len(r.seq) for r in all_input_seqs}
    seq_len_df = pd.DataFrame(seq_len_map.items(), columns=["RdRp_id", "seq_len_aa"])

    # contig/source for all IDs (tier3 not in motif_order, compute directly)
    id_meta = pd.DataFrame({
        "RdRp_id": all_ids,
        "contig":  [get_contig_name(x) for x in all_ids],
        "source":  ["lucaprot" if is_lucaprot(x) else "rdrpcatch" for x in all_ids],
    })

    # motif_order columns (tier3 rows will be NaN after left merge)
    mo_cols = motif_order[["RdRp_id", "motif_order", "n_motifs", "canonical", "permuted", "incomplete"]]

    # Combine
    summary_df = (
        tier_ser
        .merge(id_meta,        on="RdRp_id", how="left")
        .merge(seq_len_df,     on="RdRp_id", how="left")
        .merge(mo_cols,        on="RdRp_id", how="left")
        .merge(wide,           on="RdRp_id", how="left")
    )

    # Inter-motif spacings (physical order on protein)
    for ma, mb in [("A","B"), ("B","C"), ("C","D"), ("A","C"), ("A","D"), ("B","D")]:
        c1, c2 = f"mot{ma}_q1", f"mot{mb}_q1"
        if c1 in summary_df.columns and c2 in summary_df.columns:
            summary_df[f"spacing_{ma}{mb}_aa"] = (summary_df[c2] - summary_df[c1]).where(
                summary_df[c1].notna() & summary_df[c2].notna()
            )

    summary_df = summary_df.sort_values(["tier", "RdRp_id"])
    out_summary = os.path.join(args.outdir, "tier_summary.tsv")
    summary_df.to_csv(out_summary, sep="\t", index=False)
    print(f"  Written: {out_summary}  ({len(summary_df)} sequences)")

    # ── Summary ──
    print("\n=== Summary ===")
    print(f"  Sequences with ≥1 motif hit : {len(motif_order)}")
    print(f"  Complete   (all 4 motifs)   : {len(motif_order_complete)}")
    print(f"  Tier 1a  ABCD canonical     : {n_tier1a}")
    print(f"  Tier 1b  CABD depermuted    : {n_tier1b}")
    print(f"  Tier 2   3-motif partial    : {n_tier2}")
    print(f"  Tier 3   no-motif           : {n_tier3}")
    print(f"  Suspicious (D not last etc) : {len(motif_order_suspicious)}")
    print(f"  Incomplete (<4 motifs)      : {motif_order['incomplete'].sum()}")
    print(f"\nOutput files:")
    print(f"  {out_hits}")
    print(f"  {out_faa}")
    print(f"  {out_order}")
    print(f"  {out_order_complete}")
    print(f"  {out_depermed}")
    print(f"  {out_tier1a}")
    print(f"  {out_tier1b}")
    print(f"  {out_tier2}")
    print(f"  {out_tier3}")
    print(f"  {out_suspicious}")
    print(f"  {out_summary}")


if __name__ == "__main__":
    main()
