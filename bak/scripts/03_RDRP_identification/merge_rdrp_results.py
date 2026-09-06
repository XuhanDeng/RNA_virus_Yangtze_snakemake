#!/usr/bin/env python3
"""
Merge per-sample RdRp results from RdRpCATCH+palmscan and LucaProt+palmscan.

Inputs
------
--rdrpcatch         RdRpCATCH annotated TSV        (one row per contig)
--palmscan          palm_annot TSV on RC AA output  (one row per RC hit)
--lucaprot          LucaProt filtered CSV           (one row per ORF)
--lucaprot-palmscan palm_annot TSV on LP proteins   (one row per LP ORF)

Logic
-----
1. Keep only regions where palmscan confirms the palmprint motif exactly:
   pssm_ABC == "ABC" or "CAB" (case-insensitive). This applies independently
   to RC regions and LP ORFs.
2. Dedup by overlap: for each contig, any motif-confirmed LP ORF that overlaps
   (same strand, >=50% overlap of RC region length) a motif-confirmed RC
   region is dropped -- RdRpCATCH wins. Non-overlapping LP ORFs survive as
   independent candidate regions. Dropped-by-overlap LP hits are NOT
   discarded from the record -- they are tracked per contig via
   lucaprot_also_identified / lucaprot_dropped_reason, so identification
   ability can still be evaluated even when RdRpCATCH's region wins.
3. Per contig, keep only the single longest surviving region (by AA length).
   That is the contig's final RdRp call.

Output
------
--output            one row per contig with a motif-confirmed RdRp region,
                     plus lucaprot_also_identified columns recording whether
                     LucaProt independently hit this contig even if its
                     region was dropped by overlap dedup.
--output-regions    all surviving candidate regions before longest-only
                     selection (for inspection)
"""

import argparse
import os
import re
import sys

import pandas as pd


# ── constants ─────────────────────────────────────────────────────────────────

OVERLAP_THRESHOLD = 0.5
VALID_MOTIFS = {"ABC", "CAB"}

REGION_COLUMNS = [
    "contig", "source", "pssm_ABC", "aa_start", "aa_end", "aa_length",
    "frame", "strand", "seq_len_aa", "orf_start", "orf_end",
    "lucaprot_protein_id", "lucaprot_prob",
]

# extra columns recording LucaProt's independent identification of a contig,
# even when its region was dropped by overlap dedup (RdRpCATCH wins ties).
# Only present on the final per-contig table, not the region table.
LUCAPROT_EVIDENCE_COLUMNS = [
    "lucaprot_also_identified",   # True if LucaProt independently produced a
                                   # motif-confirmed ORF on this contig (kept or dropped)
    "lucaprot_overlap_dropped",   # True if that LP ORF was dropped by overlap dedup
    "lucaprot_dropped_pssm_ABC",
    "lucaprot_dropped_aa_length",
    "lucaprot_dropped_protein_id",
    "lucaprot_dropped_prob",
]


# ── helpers ───────────────────────────────────────────────────────────────────

def _motif_confirmed(val) -> bool:
    v = str(val).strip().upper() if pd.notna(val) else ""
    return v in VALID_MOTIFS


def _parse_frame(frame_str: str):
    m = re.search(r'frame=(-?\d+)', str(frame_str))
    return int(m.group(1)) if m else None


def _rc_strand(frame: int) -> str:
    return "+" if frame > 0 else "-"


def _lp_strand(orf_start: int, orf_end: int) -> str:
    return "+" if orf_start <= orf_end else "-"


def _rc_to_nt(rc_aa_start: int, rc_aa_end: int, frame: int, seq_len_aa: int):
    """Convert RC AA region to contig nucleotide coordinates."""
    fa = abs(frame)
    if frame > 0:
        nt_s = frame + (rc_aa_start - 1) * 3
        nt_e = frame + (rc_aa_end   - 1) * 3 + 2
    else:
        contig_len = seq_len_aa * 3 + fa
        nt_s = contig_len - (fa + (rc_aa_end   - 1) * 3 + 2)
        nt_e = contig_len - (fa + (rc_aa_start - 1) * 3)
    if nt_e < nt_s:
        nt_s, nt_e = nt_e, nt_s
    return nt_s, nt_e


def _overlap_fraction(rc_nt_s: int, rc_nt_e: int, lp_s: int, lp_e: int) -> float:
    lp_lo = min(lp_s, lp_e)
    lp_hi = max(lp_s, lp_e)
    ovl   = max(0, min(rc_nt_e, lp_hi) - max(rc_nt_s, lp_lo) + 1)
    rc_len = rc_nt_e - rc_nt_s + 1
    return ovl / rc_len if rc_len > 0 else 0.0


# ── loaders ───────────────────────────────────────────────────────────────────

def load_rc_regions(rdrpcatch_path: str, palmscan_path: str) -> pd.DataFrame:
    """Return motif-confirmed RC regions: one row per contig (RC gives at most one region)."""
    if not rdrpcatch_path or not os.path.exists(rdrpcatch_path) or os.path.getsize(rdrpcatch_path) == 0:
        return pd.DataFrame()
    rc = pd.read_csv(rdrpcatch_path, sep="\t", dtype=str)
    if rc.empty or "Sequence_name" not in rc.columns:
        return pd.DataFrame()
    rc = rc.rename(columns={"Sequence_name": "contig"})

    if not palmscan_path or not os.path.exists(palmscan_path) or os.path.getsize(palmscan_path) == 0:
        ps = pd.DataFrame(columns=["contig", "pssm_ABC"])
    else:
        ps = pd.read_csv(palmscan_path, sep="\t", dtype=str)
        if ps.empty or "Label" not in ps.columns:
            ps = pd.DataFrame(columns=["contig", "pssm_ABC"])
        else:
            # Label example: AnQing_D_33_0000000002_frame=-1_RdRp_1410-1766
            ps["contig"] = ps["Label"].apply(
                lambda s: re.sub(r'_frame=[^\s_]+$', '',
                          re.sub(r'_RdRp_\d+-\d+$', '', str(s))))

    # best palmscan hit per contig (highest pssm_score) if multiple
    ps_by_contig = {}
    for _, row in ps.iterrows():
        c = row["contig"]
        if c not in ps_by_contig:
            ps_by_contig[c] = row
        else:
            try:
                if float(row.get("pssm_score", 0) or 0) > \
                   float(ps_by_contig[c].get("pssm_score", 0) or 0):
                    ps_by_contig[c] = row
            except (TypeError, ValueError):
                pass

    frame_col  = "Translated_sequence_name (frame)"
    seqlen_col = "Sequence_length(AA)"

    rows = []
    for _, rc_row in rc.iterrows():
        contig = rc_row["contig"]
        ps_row = ps_by_contig.get(contig)
        pssm_abc = ps_row.get("pssm_ABC") if ps_row is not None else None
        if not _motif_confirmed(pssm_abc):
            continue

        try:
            aa_start = int(rc_row.get("RdRp_from(AA)"))
            aa_end   = int(rc_row.get("RdRp_to(AA)"))
        except (TypeError, ValueError):
            continue

        frame = _parse_frame(rc_row.get(frame_col, ""))
        try:
            seq_len_aa = int(rc_row.get(seqlen_col))
        except (TypeError, ValueError):
            seq_len_aa = None

        rows.append({
            "contig":        contig,
            "source":        "RdRpCATCH",
            "pssm_ABC":      pssm_abc,
            "aa_start":      aa_start,
            "aa_end":        aa_end,
            "aa_length":     aa_end - aa_start + 1,
            "frame":         frame,
            "strand":        _rc_strand(frame) if frame is not None else None,
            "seq_len_aa":    seq_len_aa,
            "lucaprot_protein_id": None,
            "lucaprot_prob":       None,
        })

    return pd.DataFrame(rows)


def load_lp_regions(lucaprot_path: str, lucaprot_palmscan_path: str) -> pd.DataFrame:
    """Return motif-confirmed LucaProt ORF regions: possibly multiple per contig."""
    if not lucaprot_path or not os.path.exists(lucaprot_path) or os.path.getsize(lucaprot_path) == 0:
        return pd.DataFrame()
    lp = pd.read_csv(lucaprot_path, dtype=str)
    if lp.empty or "protein_id" not in lp.columns:
        return pd.DataFrame()

    if (not lucaprot_palmscan_path or not os.path.exists(lucaprot_palmscan_path)
            or os.path.getsize(lucaprot_palmscan_path) == 0):
        lp_palm_map = {}
    else:
        ps = pd.read_csv(lucaprot_palmscan_path, sep="\t", dtype=str)
        if ps.empty or "Label" not in ps.columns or "pssm_ABC" not in ps.columns:
            lp_palm_map = {}
        else:
            lp_palm_map = dict(zip(ps["Label"].str.strip(), ps["pssm_ABC"]))

    def _parse(pid):
        s = str(pid).strip().lstrip(">")
        s = re.sub(r"^lcl\|", "", s)
        label = s.split()[0]              # ORF9_AnQing_D_33_0000000002:8171:2670
        m = re.search(r'^(.+):(\d+):(\d+)$', label)
        if not m:
            return None, label, None, None
        full_id   = m.group(1)
        orf_start = int(m.group(2))
        orf_end   = int(m.group(3))
        contig    = re.sub(r'^ORF\d+_', '', full_id)
        return contig, label, orf_start, orf_end

    rows = []
    for _, row in lp.iterrows():
        contig, label, orf_start, orf_end = _parse(row["protein_id"])
        if contig is None:
            continue
        pssm_abc = lp_palm_map.get(label)
        if not _motif_confirmed(pssm_abc):
            continue

        nt_length = abs(orf_end - orf_start) + 1
        rows.append({
            "contig":        contig,
            "source":        "LucaProt",
            "pssm_ABC":      pssm_abc,
            "orf_start":     orf_start,
            "orf_end":       orf_end,
            "aa_length":     nt_length / 3.0,
            "strand":        _lp_strand(orf_start, orf_end),
            "lucaprot_protein_id": row["protein_id"],
            "lucaprot_prob":       row.get("prob"),
        })

    return pd.DataFrame(rows)


# ── dedup by overlap (RC wins) ─────────────────────────────────────────────────

def dedup_lp_against_rc(rc_regions: pd.DataFrame, lp_regions: pd.DataFrame):
    """Drop LP regions that overlap >=50% (same strand) with an RC region on the same contig.

    Returns (kept_lp_regions, dropped_lp_regions) -- dropped_lp_regions is
    retained (not discarded) so callers can still report that LucaProt
    independently identified the contig, even though its region lost the
    overlap tie-break to RdRpCATCH.
    """
    if lp_regions.empty:
        return lp_regions, lp_regions

    rc_by_contig = {}
    for _, row in rc_regions.iterrows():
        rc_by_contig.setdefault(row["contig"], []).append(row)

    keep_mask = []
    for _, lp_row in lp_regions.iterrows():
        contig = lp_row["contig"]
        rc_candidates = rc_by_contig.get(contig, [])
        dropped = False
        for rc_row in rc_candidates:
            if rc_row.get("frame") is None or rc_row.get("seq_len_aa") is None:
                continue
            if rc_row["strand"] != lp_row["strand"]:
                continue
            rc_nt_s, rc_nt_e = _rc_to_nt(
                rc_row["aa_start"], rc_row["aa_end"], rc_row["frame"], rc_row["seq_len_aa"])
            frac = _overlap_fraction(rc_nt_s, rc_nt_e, lp_row["orf_start"], lp_row["orf_end"])
            if frac >= OVERLAP_THRESHOLD:
                dropped = True
                break
        keep_mask.append(not dropped)

    keep_mask = pd.Series(keep_mask, index=lp_regions.index)
    kept    = lp_regions[keep_mask].reset_index(drop=True)
    dropped = lp_regions[~keep_mask].reset_index(drop=True)
    return kept, dropped


# ── keep longest region per contig ────────────────────────────────────────────

def keep_longest_per_contig(regions: pd.DataFrame) -> pd.DataFrame:
    if regions.empty:
        return regions
    idx = regions.groupby("contig")["aa_length"].idxmax()
    return regions.loc[idx].reset_index(drop=True)


# ── entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Merge motif-confirmed RdRp regions (RdRpCATCH + LucaProt, palmscan-verified)."
    )
    parser.add_argument("--rdrpcatch",         required=True)
    parser.add_argument("--palmscan",          required=False, default=None)
    parser.add_argument("--lucaprot",          required=False, default=None)
    parser.add_argument("--lucaprot-palmscan", required=False, default=None)
    parser.add_argument("--output",            required=True,
                        help="Per-contig final RdRp call (longest motif-confirmed region)")
    parser.add_argument("--output-regions",    required=False, default=None,
                        help="All surviving candidate regions before longest-only selection")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    rc_regions = load_rc_regions(args.rdrpcatch, args.palmscan)
    lp_regions = load_lp_regions(args.lucaprot, args.lucaprot_palmscan)

    print(f"RC motif-confirmed regions:  {len(rc_regions)}", file=sys.stderr)
    print(f"LP motif-confirmed ORFs:     {len(lp_regions)}", file=sys.stderr)

    lp_regions_deduped, lp_regions_dropped = dedup_lp_against_rc(rc_regions, lp_regions)
    n_dropped = len(lp_regions_dropped)
    print(f"LP ORFs dropped (overlap w/ RC): {n_dropped}", file=sys.stderr)

    all_regions = pd.concat([rc_regions, lp_regions_deduped], ignore_index=True, sort=False)
    # ensure a stable column schema even when no regions survive (empty samples),
    # so downstream pd.read_csv never hits EmptyDataError on a columnless file.
    all_regions = all_regions.reindex(columns=REGION_COLUMNS)
    print(f"Surviving candidate regions: {len(all_regions)}", file=sys.stderr)

    if args.output_regions:
        os.makedirs(os.path.dirname(os.path.abspath(args.output_regions)), exist_ok=True)
        all_regions.to_csv(args.output_regions, sep="\t", index=False)
        print(f"Region table -> {args.output_regions}", file=sys.stderr)

    final = keep_longest_per_contig(all_regions)
    final = final.reindex(columns=REGION_COLUMNS)

    # attach LucaProt identification evidence, independent of which region won.
    # A contig may have LucaProt evidence via two paths:
    #   (a) its final call IS the LucaProt region (source == "LucaProt")
    #   (b) LucaProt hit it too, but that region was dropped by overlap dedup
    #       (RdRpCATCH's region won and became the final call instead)
    final["lucaprot_also_identified"] = final["source"] == "LucaProt"
    final["lucaprot_overlap_dropped"] = False
    for col in ["lucaprot_dropped_pssm_ABC", "lucaprot_dropped_aa_length",
                "lucaprot_dropped_protein_id", "lucaprot_dropped_prob"]:
        final[col] = None

    if not lp_regions_dropped.empty:
        # one dropped LP region can exist per contig per overlap check; if a
        # contig has multiple dropped LP ORFs, keep the longest for reporting
        dropped_best = (lp_regions_dropped
                         .sort_values("aa_length", ascending=False)
                         .drop_duplicates("contig", keep="first")
                         .set_index("contig"))
        for contig, drow in dropped_best.iterrows():
            mask = final["contig"] == contig
            if not mask.any():
                continue
            final.loc[mask, "lucaprot_also_identified"] = True
            final.loc[mask, "lucaprot_overlap_dropped"] = True
            final.loc[mask, "lucaprot_dropped_pssm_ABC"]    = drow.get("pssm_ABC")
            final.loc[mask, "lucaprot_dropped_aa_length"]   = drow.get("aa_length")
            final.loc[mask, "lucaprot_dropped_protein_id"]  = drow.get("lucaprot_protein_id")
            final.loc[mask, "lucaprot_dropped_prob"]        = drow.get("lucaprot_prob")

    final.to_csv(args.output, sep="\t", index=False)
    print(f"Final contig table: {len(final)} contigs -> {args.output}", file=sys.stderr)
    if not final.empty:
        print(f"  by source: {final['source'].value_counts().to_dict()}", file=sys.stderr)
        print(f"  also identified by LucaProt (kept or dropped): "
              f"{int(final['lucaprot_also_identified'].sum())}", file=sys.stderr)


if __name__ == "__main__":
    main()
