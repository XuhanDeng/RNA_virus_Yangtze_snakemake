#!/usr/bin/env python3
"""
Annotate diamond blastp hits with taxonomy and filter to Riboviria.

Since nr.dmnd was built without --taxonmap/--taxonnodes, diamond cannot filter
by --taxonlist at search time. This script does the equivalent post-hoc:
  1. Load prot.accession2taxid.tsv  (accession.version -> taxid)
  2. Load nodes.dmp                 (taxid -> parent, rank)
  3. For each diamond hit, look up the taxid of the hit accession and check
     whether Riboviria (default taxid 2559587) is in its lineage.
  4. Write a filtered TSV (same format as input + staxids column) containing
     only hits that are within Riboviria.

Input diamond TSV columns: qseqid sseqid evalue bitscore
Output TSV columns:        qseqid sseqid evalue bitscore staxids

Usage:
    python annotate_nr_taxonomy.py \
        --diamond        hits.tsv \
        --acc2tax        prot.accession2taxid.tsv \
        --nodes          taxdump/nodes.dmp \
        --riboviria-taxid 2559587 \
        --output         hits.annotated.tsv
"""

import argparse
import sys

import pandas as pd


def load_nodes(nodes_path):
    parent = {}
    with open(nodes_path) as fh:
        for line in fh:
            parts = line.split("\t|\t")
            taxid = int(parts[0].strip())
            par   = int(parts[1].strip())
            parent[taxid] = par
    return parent


def is_under(taxid, target, parent):
    visited = set()
    cur = taxid
    while cur not in visited:
        if cur == target:
            return True
        if cur == 1:
            break
        visited.add(cur)
        cur = parent.get(cur, 1)
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--diamond",          required=True)
    parser.add_argument("--acc2tax",          required=True)
    parser.add_argument("--nodes",            required=True)
    parser.add_argument("--riboviria-taxid",  type=int, default=2559587)
    parser.add_argument("--output",           required=True)
    args = parser.parse_args()

    print("Loading diamond hits...", file=sys.stderr)
    cols = ["qseqid", "sseqid", "evalue", "bitscore"]
    df = pd.read_csv(args.diamond, sep="\t", header=None, names=cols, dtype=str)
    if df.empty:
        df["staxids"] = []
        df.to_csv(args.output, sep="\t", index=False, header=False)
        print("No hits — empty output.", file=sys.stderr)
        return

    print(f"  {len(df)} hits for {df['qseqid'].nunique()} queries.", file=sys.stderr)

    # keep best hit per query (evalue ascending)
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df = df.sort_values(["qseqid", "evalue"]).drop_duplicates("qseqid", keep="first")

    # collect unique accessions needed
    accessions = set(df["sseqid"].str.split("|").str[-1])  # handle xxx|acc format
    df["_acc"] = df["sseqid"].str.split("|").str[-1]

    print("Loading accession->taxid map (this may take a few minutes)...", file=sys.stderr)
    acc2tax = {}
    with open(args.acc2tax) as fh:
        next(fh)  # skip header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0] in accessions:
                try:
                    acc2tax[parts[0]] = int(parts[1])
                except ValueError:
                    pass
    print(f"  Matched {len(acc2tax)} / {len(accessions)} accessions.", file=sys.stderr)

    print("Loading taxonomy nodes...", file=sys.stderr)
    parent = load_nodes(args.nodes)

    print("Filtering to Riboviria lineage...", file=sys.stderr)
    rows = []
    for _, row in df.iterrows():
        taxid = acc2tax.get(row["_acc"])
        if taxid is None:
            continue
        if is_under(taxid, args.riboviria_taxid, parent):
            rows.append({
                "qseqid":  row["qseqid"],
                "sseqid":  row["sseqid"],
                "evalue":  row["evalue"],
                "bitscore": row["bitscore"],
                "staxids": taxid,
            })

    out = pd.DataFrame(rows, columns=["qseqid", "sseqid", "evalue", "bitscore", "staxids"])
    out.to_csv(args.output, sep="\t", index=False, header=False)
    print(f"Kept {len(out)} / {len(df)} hits within Riboviria.", file=sys.stderr)
    print(f"  -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
