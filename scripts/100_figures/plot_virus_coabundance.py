#!/usr/bin/env python3
"""
One figure per seed (known) virus showing TPM profiles of itself and all
correlated partners across 48 samples.

Line styles:
  - Seed virus:          black, solid, linewidth=2.5
  - Known–known targets: solid, colored
  - Known–unknown targets (contigs): dashed, colored

Usage:
    python plot_virus_coabundance.py \
        --tpm       rna_virus_recalc_tpm.tsv \
        --kk        known_known_pair/r0.8.tsv \
        --ku        known_unknown_pair/r0.8.tsv \
        --outdir    result/100_figure/7_coabundance
"""

import argparse
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["svg.fonttype"] = "none"
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd


_PALETTE = [
    "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4",
    "#42D4F4", "#F032E6", "#9A6324", "#469990", "#800000",
    "#AAFFC3", "#808000", "#000075", "#FF4500", "#00CED1",
    "#FF1493", "#7FFF00", "#DC143C", "#00BFFF", "#FF8C00",
    "#9400D3", "#00FA9A", "#FF6347", "#4682B4", "#DA70D6",
]


def sample_sort_key(col: str) -> int:
    m = re.search(r"_(\d+)_tpm$", col)
    return int(m.group(1)) if m else 0


def short_label(name: str) -> str:
    """Strip t__ prefix for display."""
    return name.removeprefix("t__")


def load_tpm(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    tpm_cols = [c for c in df.columns if c.endswith("_tpm")]
    tpm_cols = sorted(tpm_cols, key=sample_sort_key)
    df = df.set_index("seq_id")[tpm_cols]
    df.columns = [c.replace("_tpm", "") for c in df.columns]
    return df


def load_pairs(pair_dir: str, threshold: float) -> dict:
    """
    Load all r{x}.tsv files in pair_dir where x >= threshold.
    Adds reverse direction (target -> seed) so pairs are mutual.
    Returns dict: seed -> list of (target, r), deduplicated (highest r kept).
    """
    import glob, re as _re
    best = {}  # (seed, target) -> r
    for path in sorted(glob.glob(os.path.join(pair_dir, "r*.tsv"))):
        m = _re.search(r"r([\d.]+)\.tsv$", path)
        if not m:
            continue
        file_r = float(m.group(1))
        if file_r < threshold:
            continue
        df = pd.read_csv(path, sep="\t")
        for _, row in df.iterrows():
            r = row["r"]
            for key in [(row["seed"], row["target"]),
                        (row["target"], row["seed"])]:
                if key not in best or r > best[key]:
                    best[key] = r
    pairs = {}
    for (seed, target), r in best.items():
        pairs.setdefault(seed, []).append((target, r))
    return pairs


def _get_series(vid: str, tpm: pd.DataFrame):
    """Return TPM series or None if missing."""
    return tpm.loc[vid] if vid in tpm.index else None


def _split_axes(lines: list, fold_threshold: float = 5.0):
    """
    Split lines into left-axis and right-axis groups.

    A line goes to the right axis if its max TPM is > fold_threshold × the
    median max (high outlier) OR < median max / fold_threshold (low outlier).
    Returns (left, right).
    """
    if not lines:
        return lines, []
    import statistics
    maxes = [s.max() for s, _, _, _, _ in lines]
    med = statistics.median(maxes)
    left, right = [], []
    for item in lines:
        mx = item[0].max()
        if med > 0 and (mx > fold_threshold * med or mx < med / fold_threshold):
            right.append(item)
        else:
            left.append(item)
    # If everything ended up on right (all similar scale), keep all on left
    if not left:
        return right, []
    return left, right


def plot_seed(seed: str, kk_targets: list, ku_targets: list,
              tpm: pd.DataFrame, outdir: str, fold_threshold: float = 5.0):
    samples = list(tpm.columns)
    x = list(range(len(samples)))

    # Build list of (series, linestyle, label, color, is_seed)
    # Colors are assigned once here so they stay consistent across left/right split
    lines = []
    color_iter = iter(_PALETTE)

    seed_s = _get_series(seed, tpm)
    if seed_s is not None:
        lines.append((seed_s, "-", f"▶ {short_label(seed)} [seed]", "black", True))
    else:
        print(f"  WARN: seed {seed} not in TPM", file=sys.stderr)

    for target, r in sorted(kk_targets, key=lambda v: -v[1]):
        s = _get_series(target, tpm)
        if s is not None:
            lines.append((s, "-",
                          f"{short_label(target)}  [kk r={r:.2f}]",
                          next(color_iter, "#888888"), False))
        else:
            print(f"  WARN: {target} not in TPM", file=sys.stderr)

    for target, r in sorted(ku_targets, key=lambda v: -v[1]):
        s = _get_series(target, tpm)
        if s is not None:
            lines.append((s, "--",
                          f"{target}  [ku r={r:.2f}]",
                          next(color_iter, "#888888"), False))
        else:
            print(f"  WARN: {target} not in TPM", file=sys.stderr)

    if not lines:
        return None

    left_lines, right_lines = _split_axes(lines, fold_threshold)

    n_partners = len(kk_targets) + len(ku_targets)
    fig_h = max(4, 3 + n_partners * 0.15)
    fig, ax = plt.subplots(figsize=(14, fig_h))
    ax2 = ax.twinx() if right_lines else None

    def draw_lines(axis, line_list, alpha=0.85):
        for series, ls, label, color, is_seed in line_list:
            lw = 2.5 if is_seed else 1.2
            zo = 5   if is_seed else 2
            axis.plot(x, series.values, color=color, linewidth=lw,
                      linestyle=ls, label=label, alpha=alpha, zorder=zo)

    draw_lines(ax, left_lines)
    if ax2:
        draw_lines(ax2, right_lines, alpha=0.7)
        ax2.set_ylabel("TPM (right axis — outlier scale)", fontsize=9, color="#666666")
        ax2.tick_params(axis="y", labelcolor="#666666", labelsize=8)
        ax2.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
        ax2.spines["top"].set_visible(False)

    ax.set_xticks(x)
    ax.set_xticklabels(samples, rotation=45, ha="right", fontsize=6)
    ax.set_ylabel("TPM", fontsize=10)
    ax.set_title(short_label(seed), fontsize=11, fontweight="bold")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.spines["top"].set_visible(False)
    if not ax2:
        ax.spines["right"].set_visible(False)

    # Combined legend — tag each label with [L] or [R] when both axes are used
    handles, labels = ax.get_legend_handles_labels()
    if ax2:
        h2, l2 = ax2.get_legend_handles_labels()
        if right_lines:
            labels = [f"{lb}  [L]" for lb in labels]
            l2     = [f"{lb}  [R]" for lb in l2]
        handles += h2
        labels  += l2

    n_total = len(handles)
    if n_total > 10:
        # Place legend below the plot so it never overlaps the lines
        fig.legend(handles, labels, fontsize=7, framealpha=0.7,
                   loc="lower center", bbox_to_anchor=(0.5, -0.02),
                   ncol=min(4, max(1, n_total // 6 + 1)),
                   bbox_transform=fig.transFigure)
    else:
        ax.legend(handles, labels, fontsize=7, loc="upper left",
                  framealpha=0.7, ncol=1)

    if right_lines:
        n_right = len(right_lines)
        note = (f"[L] left axis  |  [R] right axis — {n_right} line(s) on right "
                f"(>{fold_threshold:.0f}× or <1/{fold_threshold:.0f}× median scale)")
        fig.text(0.5, 0.01, note, ha="center", fontsize=7, color="#666666")

    plt.tight_layout()

    safe = re.sub(r'[^\w\-]', '_', seed.removeprefix("t__"))
    out_svg = os.path.join(outdir, f"{safe}.svg")
    plt.savefig(out_svg, dpi=150, bbox_inches="tight")
    plt.close()
    return out_svg


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tpm",       required=True, help="rna_virus_recalc_tpm.tsv")
    parser.add_argument("--kk-dir",    required=True, help="dir containing known_known_pair r*.tsv files")
    parser.add_argument("--ku-dir",    required=True, help="dir containing known_unknown_pair r*.tsv files")
    parser.add_argument("--threshold", type=float, default=0.7,
                        help="include all pairs with r >= threshold (default: 0.7)")
    parser.add_argument("--outdir",    required=True)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print(f"Loading TPM…", file=sys.stderr)
    tpm = load_tpm(args.tpm)

    print(f"Loading pairs (r >= {args.threshold})…", file=sys.stderr)
    kk_pairs = load_pairs(args.kk_dir, args.threshold)
    ku_pairs = load_pairs(args.ku_dir, args.threshold)

    seeds = sorted(s for s in set(list(kk_pairs.keys()) + list(ku_pairs.keys()))
                   if s.startswith("t__"))
    print(f"Seeds: {len(seeds)}", file=sys.stderr)

    for seed in seeds:
        kk_targets = kk_pairs.get(seed, [])
        ku_targets = ku_pairs.get(seed, [])
        out = plot_seed(seed, kk_targets, ku_targets, tpm, args.outdir)
        print(f"  Written → {out}  (kk={len(kk_targets)}, ku={len(ku_targets)})",
              file=sys.stderr)


if __name__ == "__main__":
    main()
