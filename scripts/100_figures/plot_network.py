#!/usr/bin/env python3
"""
Generate network figures (SVG + interactive HTML) from Spearman correlation pairs.

Two modes:
  default  : known-unknown edges filtered by --min-threshold; known-known edges same threshold
  --all-kk : known-unknown edges filtered by --min-threshold; ALL known-known edges included

Nodes  : known viruses and unknown contigs
Edges  : correlated pairs
Layout : force-directed (kamada_kawai for SVG, pyvis physics for HTML)
"""

import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams['svg.fonttype'] = 'none'   # keep text editable in Illustrator
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import pandas as pd
from pyvis.network import Network


# ── colour palette ────────────────────────────────────────────────────────────

NODE_COLORS = {
    "known":        "#4878CF",   # known reference virus
    "tier1a":       "#E84646",   # ABCD canonical — highest confidence
    "tier1b":       "#FF6B35",   # CABD depermuted
    "tier2":        "#F4C430",   # 3-motif partial
    "tier3":        "#6ACC65",   # tool-identified, <3 motifs
    "not_in_list":  "#AAAAAA",   # not an RdRp candidate
}

EDGE_COLORS = {
    "0.6": "#BDDAEF",
    "0.7": "#74B3D8",
    "0.8": "#F4A942",
    "0.9": "#E84646",
    "all": "#DDDDDD",   # known-known edges when all-kk mode
}

NODE_SIZE_SCALE  = 300
NODE_SIZE_MIN    = 100
EDGE_WIDTH_SCALE = 6.0


# ── helpers ───────────────────────────────────────────────────────────────────

def _r_group_str(val) -> str:
    try:
        return f"{float(val):.1f}"
    except (ValueError, TypeError):
        return str(val)


def _threshold_str(t: float) -> str:
    return f"{t:.1f}"


def _load_tsv(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", dtype=str)
    return df if not df.empty else pd.DataFrame()


def _update_nonzero(nonzero: dict, row, seed: str, target: str):
    try:
        nonzero[seed]   = max(nonzero.get(seed,   0),
                              int(float(row.get("seed_nonzero_samples",   0) or 0)))
        nonzero[target] = max(nonzero.get(target, 0),
                              int(float(row.get("target_nonzero_samples", 0) or 0)))
    except (ValueError, TypeError):
        pass


def build_graph(
    ku_dir: str,
    kk_dir: str,
    min_threshold: float,
    all_kk: bool,
    kk_threshold: float = None,
    kk_only: bool = False,
) -> nx.Graph:
    """
    kk_threshold : if set, load only kk pairs at exactly this threshold.
                   if None and all_kk=True, load all 4 threshold files.
                   if None and all_kk=False, load same thresholds as ku.
    """
    thresholds = [t for t in [0.6, 0.7, 0.8, 0.9] if t >= min_threshold - 1e-9]
    all_thresholds = [0.6, 0.7, 0.8, 0.9]

    G = nx.Graph()
    nonzero:  dict = {}
    rdrp_cat: dict = {}

    # ── known-unknown edges (filtered by min_threshold) ───────────────────────
    if kk_only:
        print("  kk_only mode: skipping known-unknown edges.", file=sys.stderr)
    for t in thresholds if not kk_only else []:
        fname = f"r{_threshold_str(t)}.annotated.tsv"
        df = _load_tsv(os.path.join(ku_dir, fname))
        if df.empty:
            continue
        for _, row in df.iterrows():
            seed   = str(row.get("seed",   "")).strip()
            target = str(row.get("target", "")).strip()
            if not seed or not target:
                continue
            r_val = float(row.get("r", 0))
            rg    = _r_group_str(row.get("r_group", _threshold_str(t)))
            cat   = str(row.get("rdrp_category", "not_in_list")).strip()

            if seed not in G:
                G.add_node(seed, node_type="known")
            if target not in G:
                G.add_node(target, node_type="unknown", rdrp_category=cat)
            rdrp_cat[target] = cat
            _update_nonzero(nonzero, row, seed, target)
            G.add_edge(seed, target, r=r_val, r_group=rg, edge_type="known_unknown")

    # ── known-known edges ─────────────────────────────────────────────────────
    if kk_threshold is not None:
        kk_thresholds = [kk_threshold]
        kk_label = _threshold_str(kk_threshold)
    elif all_kk:
        kk_thresholds = all_thresholds
        kk_label = "all"
    else:
        kk_thresholds = thresholds
        kk_label = None   # use per-row r_group

    for t in kk_thresholds:
        fname = f"r{_threshold_str(t)}.tsv"
        df = _load_tsv(os.path.join(kk_dir, fname))
        if df.empty:
            continue
        for _, row in df.iterrows():
            seed   = str(row.get("seed",   "")).strip()
            target = str(row.get("target", "")).strip()
            if not seed or not target:
                continue
            r_val = float(row.get("r", 0))
            rg    = kk_label if kk_label else _r_group_str(row.get("r_group", _threshold_str(t)))

            for n in (seed, target):
                if n not in G:
                    G.add_node(n, node_type="known")
            _update_nonzero(nonzero, row, seed, target)
            G.add_edge(seed, target, r=r_val, r_group=rg, edge_type="known_known")

    # attach nonzero / rdrp_category to nodes
    for n in G.nodes():
        G.nodes[n]["nonzero_samples"] = nonzero.get(n, 1)
        if G.nodes[n].get("node_type") == "unknown":
            G.nodes[n]["rdrp_category"] = rdrp_cat.get(n, "not_in_list")

    return G


def filter_not_in_list(G: nx.Graph) -> nx.Graph:
    """Remove unknown nodes with rdrp_category == 'not_in_list' and their edges."""
    remove = [n for n, d in G.nodes(data=True)
              if d.get("node_type") == "unknown"
              and d.get("rdrp_category") == "not_in_list"]
    G.remove_nodes_from(remove)
    print(f"  Removed {len(remove)} not_in_list nodes.", file=sys.stderr)
    return G


# ── node / edge styling ───────────────────────────────────────────────────────

def _node_color(G: nx.Graph, node: str) -> str:
    data = G.nodes[node]
    if data.get("node_type") == "known":
        return NODE_COLORS["known"]
    return NODE_COLORS.get(data.get("rdrp_category", "not_in_list"),
                           NODE_COLORS["not_in_list"])


def _node_size_mpl(G: nx.Graph, node: str) -> float:
    n = G.nodes[node].get("nonzero_samples", 1)
    return max(NODE_SIZE_MIN, n * NODE_SIZE_SCALE / 10)


def _edge_color(G: nx.Graph, u: str, v: str) -> str:
    rg = str(G.edges[u, v].get("r_group", "0.6"))
    return EDGE_COLORS.get(rg, "#CCCCCC")


def _edge_width(G: nx.Graph, u: str, v: str) -> float:
    r = float(G.edges[u, v].get("r", 0.6))
    return max(0.5, (r - 0.5) / 0.5 * EDGE_WIDTH_SCALE)


# ── static SVG ────────────────────────────────────────────────────────────────

def draw_svg(G: nx.Graph, out_path: str, with_label: bool,
             min_threshold: float, all_kk: bool):
    if len(G.nodes) == 0:
        print(f"  Graph empty — skipping {out_path}", file=sys.stderr)
        return

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    print(f"  Layout ({len(G.nodes)} nodes, {len(G.edges)} edges)…",
          file=sys.stderr, flush=True)

    if len(G.nodes) <= 500:
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.spring_layout(G, seed=42, k=1.5 / (len(G.nodes) ** 0.5))

    fig, ax = plt.subplots(figsize=(18, 14))
    ax.set_axis_off()

    node_list   = list(G.nodes())
    node_colors = [_node_color(G, n) for n in node_list]
    node_sizes  = [_node_size_mpl(G, n) for n in node_list]

    edge_list   = list(G.edges())
    edge_colors = [_edge_color(G, u, v) for u, v in edge_list]
    edge_widths = [_edge_width(G, u, v) for u, v in edge_list]

    nx.draw_networkx_edges(G, pos, edgelist=edge_list,
                           edge_color=edge_colors, width=edge_widths,
                           alpha=0.7, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=node_list,
                           node_color=node_colors, node_size=node_sizes,
                           alpha=0.9, ax=ax)
    if with_label:
        nx.draw_networkx_labels(G, pos, font_size=5, ax=ax)

    node_patches = [mpatches.Patch(color=v, label=k) for k, v in NODE_COLORS.items()]
    edge_entries = (
        [mpatches.Patch(color=EDGE_COLORS["all"], label="known-known (all r)")]
        if all_kk else
        [mpatches.Patch(color=v, label=f"r ≥ {k}") for k, v in EDGE_COLORS.items()
         if k != "all"]
    )
    ax.legend(handles=node_patches + edge_entries,
              loc="upper left", fontsize=7, framealpha=0.8)

    title = f"Spearman network  (r ≥ {min_threshold}"
    title += ", all known-known)" if all_kk else ")"
    ax.set_title(title, fontsize=12)

    fig.tight_layout()
    fig.savefig(out_path, format="svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved SVG: {out_path}", file=sys.stderr)


# ── interactive HTML ──────────────────────────────────────────────────────────

def draw_html(G: nx.Graph, out_path: str, with_label: bool,
              min_threshold: float, all_kk: bool):
    if len(G.nodes) == 0:
        print(f"  Graph empty — skipping {out_path}", file=sys.stderr)
        return

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)

    net = Network(height="900px", width="100%", bgcolor="#1a1a2e",
                  font_color="white", notebook=False)
    net.barnes_hut(gravity=-8000, central_gravity=0.3,
                   spring_length=200, spring_strength=0.05, damping=0.9)

    for node in G.nodes():
        data  = G.nodes[node]
        color = _node_color(G, node)
        size  = max(10, data.get("nonzero_samples", 1) * 2)
        label = node if with_label else ""
        title = (
            f"<b>{node}</b><br>"
            f"Type: {data.get('node_type', '?')}<br>"
            f"RdRP category: {data.get('rdrp_category', 'known')}<br>"
            f"Non-zero samples: {data.get('nonzero_samples', '?')}"
        )
        net.add_node(node, label=label, color=color, size=size, title=title)

    for u, v, edata in G.edges(data=True):
        color = _edge_color(G, u, v)
        width = _edge_width(G, u, v)
        title = f"r = {float(edata.get('r', 0)):.4f}  (r_group: {edata.get('r_group', '?')})"
        net.add_edge(u, v, color=color, width=width, title=title)

    net.set_options("""
    {
      "physics": {
        "enabled": true,
        "stabilization": { "iterations": 200 }
      }
    }
    """)

    net.save_graph(out_path)
    print(f"  Saved HTML: {out_path}", file=sys.stderr)


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Plot Spearman correlation network (SVG + HTML)."
    )
    parser.add_argument("--ku-dir",        required=True)
    parser.add_argument("--kk-dir",        required=True)
    parser.add_argument("--min-threshold", required=True, type=float)
    parser.add_argument("--out-svg-label",    required=True)
    parser.add_argument("--out-svg-nolabel",  required=True)
    parser.add_argument("--out-html-label",   required=True)
    parser.add_argument("--out-html-nolabel", required=True)
    parser.add_argument("--all-kk", action="store_true",
                        help="Include all known-known pairs regardless of r threshold")
    parser.add_argument("--kk-threshold", type=float, default=None,
                        help="Load known-known pairs at exactly this threshold only "
                             "(e.g. 0.6). Overrides --all-kk if both given.")
    parser.add_argument("--hide-not-in-list", action="store_true",
                        help="Remove unknown nodes with rdrp_category=not_in_list")
    parser.add_argument("--kk-only", action="store_true",
                        help="Skip all known-unknown edges; show only known-known pairs")
    args = parser.parse_args()

    kk_thr = args.kk_threshold
    print(f"Building graph (r >= {args.min_threshold}, all_kk={args.all_kk}, "
          f"kk_threshold={kk_thr}, hide_not_in_list={args.hide_not_in_list})…",
          file=sys.stderr, flush=True)
    G = build_graph(args.ku_dir, args.kk_dir, args.min_threshold,
                    args.all_kk, kk_threshold=kk_thr, kk_only=args.kk_only)

    if args.hide_not_in_list:
        G = filter_not_in_list(G)
    print(f"  Nodes: {len(G.nodes)}  Edges: {len(G.edges)}", file=sys.stderr)

    draw_svg( G, args.out_svg_label,    with_label=True,  min_threshold=args.min_threshold, all_kk=args.all_kk)
    draw_svg( G, args.out_svg_nolabel,  with_label=False, min_threshold=args.min_threshold, all_kk=args.all_kk)
    draw_html(G, args.out_html_label,   with_label=True,  min_threshold=args.min_threshold, all_kk=args.all_kk)
    draw_html(G, args.out_html_nolabel, with_label=False, min_threshold=args.min_threshold, all_kk=args.all_kk)


if __name__ == "__main__":
    main()
