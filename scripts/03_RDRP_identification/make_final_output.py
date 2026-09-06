#!/usr/bin/env python3
"""
Collect final RdRP outputs into two subfolders:

  10_final/all/          all sequences before NR filter
  10_final/nr_filtered/  NR-confirmed sequences only

Also writes final_rdrp_merged.tsv and final_contigs.fasta at the top level,
filtered to NR-confirmed contigs only.
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

import pandas as pd


def add_motif_arrangement(df: pd.DataFrame) -> pd.DataFrame:
    """Add motif_arrangement column (ABC, CAB, or other) based on motif start positions."""
    def arrangement(row):
        try:
            a = int(row["motif_a_start"])
            b = int(row["motif_b_start"])
            c = int(row["motif_c_start"])
        except (ValueError, TypeError):
            return "unknown"
        if a < b < c:
            return "ABC"
        if c < a < b:
            return "CAB"
        return "other"
    df = df.copy()
    df.insert(df.columns.get_loc("flag"), "motif_arrangement",
              df.apply(arrangement, axis=1))
    return df


def read_faa_ids(path: str) -> set:
    ids = set()
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.add(line[1:].split()[0])
    return ids


def combine_fasta(src_a: str, src_b: str, dst: str):
    with open(dst, "w") as fout:
        for src in (src_a, src_b):
            with open(src) as fin:
                fout.write(fin.read())


def combine_tsv(src_a: str, label_a: str, src_b: str, label_b: str, dst: str):
    df_a = pd.read_csv(src_a, sep="\t", dtype=str)
    df_b = pd.read_csv(src_b, sep="\t", dtype=str)
    if "source" not in df_a.columns:
        df_a.insert(0, "source", label_a)
    else:
        df_a["source"] = label_a
    if "source" not in df_b.columns:
        df_b.insert(0, "source", label_b)
    else:
        df_b["source"] = label_b
    pd.concat([df_a, df_b], ignore_index=True).to_csv(dst, sep="\t", index=False)


def filter_fasta(src: str, keep_ids: set, dst: str):
    kept = 0
    with open(src) as fin, open(dst, "w") as fout:
        write = False
        for line in fin:
            if line.startswith(">"):
                write = line[1:].split()[0] in keep_ids
                if write:
                    kept += 1
            if write:
                fout.write(line)
    return kept


def to_contig(pid: str) -> str:
    pid = re.sub(r"^ORF\d+_", "", pid)
    pid = pid.split(":")[0]
    pid = re.sub(r"_frame=.*", "", pid)
    return pid


def main():
    p = argparse.ArgumentParser(description=__doc__)
    # NR-confirmed proteins (from filter_nr_viral_origin)
    p.add_argument("--rc-confirmed",   required=True, help="RdRpCATCH_viral_confirmed.faa")
    p.add_argument("--lp-confirmed",   required=True, help="LucaProt_viral_confirmed.faa")
    # palm region files (all sequences, pre-NR filter)
    p.add_argument("--rc-full",        required=True)
    p.add_argument("--lp-full",        required=True)
    p.add_argument("--rc-core",        required=True)
    p.add_argument("--lp-core",        required=True)
    p.add_argument("--rc-extended",    required=True)
    p.add_argument("--lp-extended",    required=True)
    p.add_argument("--rc-palm-tsv",    required=True)
    p.add_argument("--lp-palm-tsv",    required=True)
    # other inputs
    p.add_argument("--merged",         required=True)
    p.add_argument("--contigs",        required=True)
    # output dirs
    p.add_argument("--outdir",         required=True)
    args = p.parse_args()

    outdir     = Path(args.outdir)
    all_dir    = outdir / "all"
    nr_dir     = outdir / "nr_filtered"
    all_dir.mkdir(parents=True, exist_ok=True)
    nr_dir.mkdir(parents=True, exist_ok=True)

    # ── all/ : copy pre-NR palm files directly ────────────────────────────────
    for src, name in [
        (args.rc_full,     "RdRpCATCH_full.faa"),
        (args.lp_full,     "LucaProt_full.faa"),
        (args.rc_core,     "RdRpCATCH_palm_core.faa"),
        (args.lp_core,     "LucaProt_palm_core.faa"),
        (args.rc_extended, "RdRpCATCH_palm_extended.faa"),
        (args.lp_extended, "LucaProt_palm_extended.faa"),
    ]:
        shutil.copy(src, all_dir / name)
    for src, name in [
        (args.rc_palm_tsv, "RdRpCATCH_palm_regions.tsv"),
        (args.lp_palm_tsv, "LucaProt_palm_regions.tsv"),
    ]:
        df = pd.read_csv(src, sep="\t", dtype=str)
        add_motif_arrangement(df).to_csv(all_dir / name, sep="\t", index=False)
    print(f"all/ written: {all_dir}", file=sys.stderr)

    # ── NR-confirmed IDs ──────────────────────────────────────────────────────
    rc_ids = read_faa_ids(args.rc_confirmed)
    lp_ids = read_faa_ids(args.lp_confirmed)
    print(f"NR-confirmed — RdRpCATCH: {len(rc_ids)}  LucaProt: {len(lp_ids)}",
          file=sys.stderr)

    # ── nr_filtered/ : confirmed faa + filtered palm files ───────────────────
    shutil.copy(args.rc_confirmed, nr_dir / "RdRpCATCH.faa")
    shutil.copy(args.lp_confirmed, nr_dir / "LucaProt.faa")

    for src, ids, name in [
        (args.rc_full,     rc_ids, "RdRpCATCH_full.faa"),
        (args.lp_full,     lp_ids, "LucaProt_full.faa"),
        (args.rc_core,     rc_ids, "RdRpCATCH_palm_core.faa"),
        (args.lp_core,     lp_ids, "LucaProt_palm_core.faa"),
        (args.rc_extended, rc_ids, "RdRpCATCH_palm_extended.faa"),
        (args.lp_extended, lp_ids, "LucaProt_palm_extended.faa"),
    ]:
        n = filter_fasta(src, ids, str(nr_dir / name))
        print(f"  {name}: {n} sequences", file=sys.stderr)

    # filter palm TSVs by confirmed IDs + add motif_arrangement
    for src, ids, name in [
        (args.rc_palm_tsv, rc_ids, "RdRpCATCH_palm_regions.tsv"),
        (args.lp_palm_tsv, lp_ids, "LucaProt_palm_regions.tsv"),
    ]:
        df = pd.read_csv(src, sep="\t", dtype=str)
        filtered = add_motif_arrangement(df[df["seq_id"].isin(ids)])
        filtered.to_csv(nr_dir / name, sep="\t", index=False)
        print(f"  {name}: {len(filtered)} rows", file=sys.stderr)

    print(f"nr_filtered/ written: {nr_dir}", file=sys.stderr)

    # ── combined files (RdRpCATCH + LucaProt merged) ──────────────────────────
    for subdir, rc_tsv, lp_tsv in [
        (all_dir, all_dir / "RdRpCATCH_palm_regions.tsv", all_dir / "LucaProt_palm_regions.tsv"),
        (nr_dir,  nr_dir  / "RdRpCATCH_palm_regions.tsv", nr_dir  / "LucaProt_palm_regions.tsv"),
    ]:
        combine_fasta(str(subdir / "RdRpCATCH_full.faa"),
                      str(subdir / "LucaProt_full.faa"),
                      str(subdir / "combined_full.faa"))
        combine_fasta(str(subdir / "RdRpCATCH_palm_core.faa"),
                      str(subdir / "LucaProt_palm_core.faa"),
                      str(subdir / "combined_palm_core.faa"))
        combine_fasta(str(subdir / "RdRpCATCH_palm_extended.faa"),
                      str(subdir / "LucaProt_palm_extended.faa"),
                      str(subdir / "combined_palm_extended.faa"))
        combine_tsv(str(rc_tsv), "RdRpCATCH",
                    str(lp_tsv), "LucaProt",
                    str(subdir / "combined_palm_regions.tsv"))
    # combined NR-confirmed proteins (the .faa equivalent of RdRpCATCH.faa + LucaProt.faa)
    combine_fasta(str(nr_dir / "RdRpCATCH.faa"),
                  str(nr_dir / "LucaProt.faa"),
                  str(nr_dir / "combined.faa"))
    print(f"combined files written in all/ and nr_filtered/", file=sys.stderr)

    # ── top-level: merged TSV + contigs filtered to NR-confirmed ─────────────
    confirmed_contigs = {to_contig(pid) for pid in rc_ids | lp_ids}
    print(f"Confirmed contigs: {len(confirmed_contigs)}", file=sys.stderr)

    df = pd.read_csv(args.merged, sep="\t", dtype=str)
    id_col = next((c for c in ("contig", "contig_id", "contigID", "seq_id", "seqid", "query")
                   if c in df.columns), df.columns[0])
    filtered = df[df[id_col].isin(confirmed_contigs)]
    filtered.to_csv(outdir / "final_rdrp_merged.tsv", sep="\t", index=False)
    print(f"final_rdrp_merged.tsv: {len(filtered)} / {len(df)} rows", file=sys.stderr)

    kept = filter_fasta(args.contigs, confirmed_contigs,
                        str(outdir / "final_contigs.fasta"))
    print(f"final_contigs.fasta: {kept} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
