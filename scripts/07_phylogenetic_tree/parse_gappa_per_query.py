#!/usr/bin/env python3
"""
Parse gappa examine assign per_query.tsv and produce:
  1. phylum_annotation.tsv  — one row per query, phylum-level (afract >= threshold)
  2. full_taxonomy_annotation.tsv — one row per query, deepest confident taxopath

Usage:
    python parse_gappa_per_query.py --input per_query.tsv --outdir <dir> [--threshold 0.66]
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def parse_per_query(path: str):
    """
    Return dict: query -> list of row dicts (one per taxopath entry).
    Rows within each query are sorted by taxopath depth (ascending).
    """
    queries = defaultdict(list)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            queries[row["name"]].append({
                "LWR":     float(row["LWR"]),
                "fract":   float(row["fract"]),
                "aLWR":    float(row["aLWR"]),
                "afract":  float(row["afract"]),
                "taxopath": row["taxopath"].strip(),
            })
    # Sort each query's rows by taxopath depth, then by taxopath string for stability
    for name in queries:
        queries[name].sort(key=lambda r: (r["taxopath"].count(";"), r["taxopath"]))
    return queries


def build_phylum_annotation(queries: dict, threshold: float) -> list:
    """File 1: one row per query at phylum level."""
    rows = []
    for query, entries in queries.items():
        # Phylum rows = depth 0 (no semicolon)
        phylum_rows = [e for e in entries if ";" not in e["taxopath"]]

        if not phylum_rows:
            rows.append({
                "query":   query,
                "phylum":  "Unassigned",
                "LWR":     "",
                "fract":   "",
                "aLWR":    "",
                "afract":  "",
            })
            continue

        # Among phylum rows, pick any with afract >= threshold
        confident = [e for e in phylum_rows if e["afract"] >= threshold]

        if confident:
            # Pick highest afract among confident phylum rows
            best = max(confident, key=lambda e: e["afract"])
            rows.append({
                "query":   query,
                "phylum":  best["taxopath"],
                "LWR":     f"{best['LWR']:.6g}",
                "fract":   f"{best['fract']:.6g}",
                "aLWR":    f"{best['aLWR']:.6g}",
                "afract":  f"{best['afract']:.6g}",
            })
        else:
            # Fallback: best phylum row by afract, mark Unassigned
            best = max(phylum_rows, key=lambda e: e["afract"])
            rows.append({
                "query":   query,
                "phylum":  "Unassigned",
                "LWR":     f"{best['LWR']:.6g}",
                "fract":   f"{best['fract']:.6g}",
                "aLWR":    f"{best['aLWR']:.6g}",
                "afract":  f"{best['afract']:.6g}",
            })

    return rows


def build_full_taxonomy_annotation(queries: dict, threshold: float) -> list:
    """File 2: one row per query, deepest taxopath where afract >= threshold."""
    rows = []
    for query, entries in queries.items():
        # Best single placement (highest LWR)
        best_lwr_row = max(entries, key=lambda e: e["LWR"])
        full_taxopath_best = best_lwr_row["taxopath"]

        # Deepest confident: walk from deepest to shallowest
        # entries are sorted by depth asc; reverse for deepest-first search
        confident_entry = None
        for e in reversed(entries):
            if e["afract"] >= threshold:
                confident_entry = e
                break

        if confident_entry is None:
            rows.append({
                "query":                query,
                "confident_taxopath":   "Unassigned",
                "deepest_rank_depth":   0,
                "afract_at_deepest":    "",
                "full_taxopath_best":   full_taxopath_best,
            })
        else:
            tp = confident_entry["taxopath"]
            depth = tp.count(";") + 1
            rows.append({
                "query":                query,
                "confident_taxopath":   tp,
                "deepest_rank_depth":   depth,
                "afract_at_deepest":    f"{confident_entry['afract']:.6g}",
                "full_taxopath_best":   full_taxopath_best,
            })

    return rows


def write_tsv(path: str, fieldnames: list, rows: list):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Written → {path}  ({len(rows)} rows)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",     required=True, help="per_query.tsv from gappa examine assign")
    parser.add_argument("--outdir",    required=True, help="output directory")
    parser.add_argument("--threshold", type=float, default=0.66,
                        help="afract threshold for confident assignment (default: 0.66)")
    args = parser.parse_args()

    queries = parse_per_query(args.input)
    print(f"Queries loaded: {len(queries)}", file=sys.stderr)

    phylum_rows = build_phylum_annotation(queries, args.threshold)
    full_rows   = build_full_taxonomy_annotation(queries, args.threshold)

    write_tsv(
        os.path.join(args.outdir, "phylum_annotation.tsv"),
        ["query", "phylum", "LWR", "fract", "aLWR", "afract"],
        phylum_rows,
    )
    write_tsv(
        os.path.join(args.outdir, "full_taxonomy_annotation.tsv"),
        ["query", "confident_taxopath", "deepest_rank_depth", "afract_at_deepest", "full_taxopath_best"],
        full_rows,
    )

    # Summary stats
    assigned_phylum = sum(1 for r in phylum_rows if r["phylum"] != "Unassigned")
    assigned_full   = sum(1 for r in full_rows   if r["confident_taxopath"] != "Unassigned")
    print(f"Phylum assigned (afract≥{args.threshold}): {assigned_phylum}/{len(phylum_rows)}", file=sys.stderr)
    print(f"Deep confident  (afract≥{args.threshold}): {assigned_full}/{len(full_rows)}", file=sys.stderr)


if __name__ == "__main__":
    main()
