#!/usr/bin/env python3
"""
Classify contig RdRp tips by their nearest ICTV reference tip in the phylogenetic tree.

Strategy:
  For each contig tip, walk up the tree until the current clade contains at least one
  ICTV tip with VMR taxonomy. Among those ICTV tips, find the one with the smallest
  patristic distance (sum of branch lengths through their lowest common ancestor).

Inputs:
  --tree        renamed Newick tree (my_tree.renamed.nwk)
  --vmr         ICTV VMR Excel file (VMR_*.xlsx)
  --lookup      my_tree.renamed_lookup.tsv (to identify contig vs ICTV tips)
  --outdir      output directory

Output:
  contig_classification.tsv  — one row per contig tip
"""

import argparse
import os
import re
import sys

try:
    from Bio import Phylo
except ImportError:
    sys.exit("ERROR: biopython not installed. Run: conda install -c conda-forge biopython")

try:
    import openpyxl
except ImportError:
    sys.exit("ERROR: openpyxl not installed.")


# ── VMR loading ───────────────────────────────────────────────────────────────

def load_vmr(vmr_path, vmr_sheet=None):
    wb = openpyxl.load_workbook(vmr_path, read_only=True)
    if vmr_sheet and vmr_sheet in wb.sheetnames:
        ws = wb[vmr_sheet]
    else:
        ws = next((wb[s] for s in wb.sheetnames if "VMR" in s or "MSL" in s), wb.active)

    rows = list(ws.iter_rows(values_only=True))
    header = list(rows[0])

    def col(name):
        return header.index(name)

    ord_col = col("Order")
    fam_col = col("Family")
    gen_col = col("Genus")
    spe_col = col("Species")
    acc_col = col("Virus GENBANK accession")

    acc2tax = {}
    for r in rows[1:]:
        raw_acc = r[acc_col]
        if not raw_acc:
            continue
        tax = {
            "Order":   r[ord_col]  or "",
            "Family":  r[fam_col]  or "",
            "Genus":   r[gen_col]  or "",
            "Species": r[spe_col]  or "",
        }
        for acc in re.findall(r"[A-Z]+\d+", str(raw_acc)):
            acc2tax[acc.split(".")[0]] = tax

    print(f"VMR: {len(acc2tax)} accessions with taxonomy", file=sys.stderr)
    return acc2tax


# ── Lookup loading ────────────────────────────────────────────────────────────

def load_lookup(lookup_path):
    ictv_tips = set()
    contig_tips = set()
    with open(lookup_path) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            renamed = parts[1]
            if renamed.startswith("ICTV_"):
                ictv_tips.add(renamed)
            elif not renamed.startswith("esvirtu_") and not renamed.startswith("rt."):
                contig_tips.add(renamed)
    print(f"Lookup: {len(ictv_tips)} ICTV, {len(contig_tips)} contigs", file=sys.stderr)
    return ictv_tips, contig_tips


def ictv_acc_base(tip_name):
    """ICTV_AB000906.1 → AB000906"""
    return tip_name[5:].split(".")[0]


# ── Fast tree traversal ───────────────────────────────────────────────────────

def build_tree_index(tree):
    """
    Pre-compute for every clade:
      - depth_from_root (sum of branch lengths from root)
      - parent map
      - terminal sets per internal node (lazily via walk)

    Returns:
      node_depth: {clade_id: float}   (clade_id = id(clade))
      parent:     {clade_id: clade}
      tip_by_name: {name: clade}
    """
    node_depth = {}
    parent = {}
    tip_by_name = {}

    def _walk(clade, depth):
        cid = id(clade)
        node_depth[cid] = depth
        if clade.is_terminal():
            tip_by_name[clade.name] = clade
        else:
            for child in clade.clades:
                parent[id(child)] = clade
                _walk(child, depth + (child.branch_length or 0.0))

    root = tree.root
    node_depth[id(root)] = 0.0
    for child in root.clades:
        parent[id(child)] = root
        _walk(child, child.branch_length or 0.0)

    return node_depth, parent, tip_by_name


def patristic(tip_a, tip_b, node_depth, parent):
    """
    Fast patristic distance via LCA:
      dist(A,B) = depth(A) + depth(B) - 2*depth(LCA)
    """
    # Collect ancestors of A (including A itself)
    ancestors_a = {}
    node = tip_a
    while node is not None:
        ancestors_a[id(node)] = node
        node = parent.get(id(node))

    # Walk up from B until we hit an ancestor of A
    node = tip_b
    while node is not None:
        if id(node) in ancestors_a:
            lca = node
            break
        node = parent.get(id(node))
    else:
        return float("inf")

    return (node_depth[id(tip_a)] + node_depth[id(tip_b)]
            - 2.0 * node_depth[id(lca)])


# ── Main classification ───────────────────────────────────────────────────────

def classify_contigs(tree, ictv_tips, contig_tips, acc2tax):
    print("Building tree index ...", file=sys.stderr)
    node_depth, parent, tip_by_name = build_tree_index(tree)

    # Map ICTV tip → taxonomy
    ictv_tax = {}
    for tip in ictv_tips:
        acc = ictv_acc_base(tip)
        ictv_tax[tip] = acc2tax.get(acc)  # None if no VMR match

    ictv_with_tax = {t for t, tax in ictv_tax.items() if tax is not None}
    print(f"ICTV tips with VMR taxonomy: {len(ictv_with_tax)}", file=sys.stderr)
    print(f"ICTV tips without VMR match: {len(ictv_tips) - len(ictv_with_tax)}", file=sys.stderr)

    # Pre-fetch clade objects for ICTV tips with tax
    ictv_clade = {t: tip_by_name[t] for t in ictv_with_tax if t in tip_by_name}

    results = []
    contigs_in_tree = contig_tips & set(tip_by_name.keys())
    print(f"Contigs in tree: {len(contigs_in_tree)} / {len(contig_tips)}", file=sys.stderr)

    for i, contig in enumerate(sorted(contigs_in_tree)):
        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(contigs_in_tree)} ...", file=sys.stderr)

        contig_clade = tip_by_name.get(contig)
        if contig_clade is None:
            continue

        # Walk up from contig; stop at the first ancestor whose terminal set
        # intersects ictv_with_tax
        candidate_ictv = []
        node = contig_clade
        while node is not None:
            terminals = {c.name for c in node.get_terminals()}
            candidates = terminals & ictv_with_tax
            if candidates:
                candidate_ictv = list(candidates)
                break
            node = parent.get(id(node))

        if not candidate_ictv:
            candidate_ictv = list(ictv_with_tax)  # whole-tree fallback

        # Find nearest among candidates
        best_tip = None
        best_dist = float("inf")
        for ictv_tip in candidate_ictv:
            ic = ictv_clade.get(ictv_tip)
            if ic is None:
                continue
            d = patristic(contig_clade, ic, node_depth, parent)
            if d < best_dist:
                best_dist = d
                best_tip = ictv_tip

        if best_tip is None:
            results.append({
                "contig": contig, "nearest_ictv": "", "ictv_accession": "",
                "Order": "", "Family": "", "Genus": "", "Species": "",
                "distance": "", "note": "no_ictv_found",
            })
            continue

        tax = ictv_tax[best_tip]
        results.append({
            "contig": contig,
            "nearest_ictv": best_tip,
            "ictv_accession": ictv_acc_base(best_tip),
            "Order":   tax["Order"],
            "Family":  tax["Family"],
            "Genus":   tax["Genus"],
            "Species": tax["Species"],
            "distance": f"{best_dist:.6f}",
            "note": "",
        })

    return results


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree",      required=True)
    parser.add_argument("--vmr",       required=True)
    parser.add_argument("--vmr-sheet", default=None)
    parser.add_argument("--lookup",    required=True)
    parser.add_argument("--outdir",    required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Loading VMR taxonomy ...", file=sys.stderr)
    acc2tax = load_vmr(args.vmr, getattr(args, "vmr_sheet", None))

    print("Loading lookup ...", file=sys.stderr)
    ictv_tips, contig_tips = load_lookup(args.lookup)

    print("Reading tree ...", file=sys.stderr)
    tree = Phylo.read(args.tree, "newick")
    print(f"Tree: {tree.count_terminals()} terminals", file=sys.stderr)

    results = classify_contigs(tree, ictv_tips, contig_tips, acc2tax)

    outpath = os.path.join(args.outdir, "contig_classification.tsv")
    cols = ["contig", "nearest_ictv", "ictv_accession",
            "Order", "Family", "Genus", "Species", "distance", "note"]
    with open(outpath, "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in results:
            out.write("\t".join(str(r[c]) for c in cols) + "\n")

    print(f"\nWritten {len(results)} rows → {outpath}", file=sys.stderr)

    from collections import Counter
    fam_counts = Counter(r["Family"] for r in results if r["Family"])
    print("\nTop families assigned:", file=sys.stderr)
    for fam, cnt in fam_counts.most_common(15):
        print(f"  {fam}: {cnt}", file=sys.stderr)


if __name__ == "__main__":
    main()
