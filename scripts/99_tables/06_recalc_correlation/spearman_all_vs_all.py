#!/usr/bin/env python3
"""
All-to-all Spearman correlation — filtered output only.
Adapted from scripts/05_esviritu/spearman_all_vs_all.py for use with
the recalculated RNA-virus-only RPKMF/TPM tables (rna_virus_recalc_rpkmf.tsv /
rna_virus_recalc_tpm.tsv).

Change from original: sample columns end with a configurable suffix
(--value-suffix, default "_rpkmf") instead of "_RPKMF".
"""

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing.shared_memory import SharedMemory
from typing import List

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as t_dist


def _parse_thresholds(text: str) -> List[float]:
    vals = [float(p.strip()) for p in text.split(",") if p.strip()]
    if not vals:
        raise SystemExit("No valid thresholds provided.")
    return sorted(set(vals))


def _label_col(df: pd.DataFrame, value_suffix: str) -> str:
    for c in df.columns:
        if not c.endswith(value_suffix):
            return c
    raise SystemExit("No label column found in input.")


def _prerank(data: np.ndarray):
    ranked = np.apply_along_axis(rankdata, 1, data).astype(float)
    ranked -= ranked.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(ranked, axis=1)
    return ranked, norms


def _spearman_block(ranked_a, norms_a, ranked_b, norms_b, m):
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (ranked_a @ ranked_b.T) / np.outer(norms_a, norms_b)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((m - 2) / (1 - r ** 2))
    p = 2 * t_dist.sf(np.abs(t), df=m - 2)
    return r, p


def _worker(
    chunk_indices: np.ndarray,
    all_indices: np.ndarray,
    shm_data_name: str, data_shape: tuple, data_dtype,
    shm_ranked_name: str, ranked_shape: tuple,
    shm_norms_name: str, norms_shape: tuple,
    m: int,
    labels_chunk: np.ndarray,
    labels_all: np.ndarray,
    min_thr: float,
    p_thr: float,
    thresholds: List[float],
    nonzero_map: dict,
):
    shm_data   = SharedMemory(name=shm_data_name,   create=False)
    shm_ranked = SharedMemory(name=shm_ranked_name, create=False)
    shm_norms  = SharedMemory(name=shm_norms_name,  create=False)

    data        = np.ndarray(data_shape,   dtype=data_dtype,  buffer=shm_data.buf)
    ranked_all  = np.ndarray(ranked_shape, dtype=np.float64,  buffer=shm_ranked.buf)
    norms_all   = np.ndarray(norms_shape,  dtype=np.float64,  buffer=shm_norms.buf)

    ranked_chunk, norms_chunk = _prerank(data[chunk_indices])
    r_block, p_block = _spearman_block(ranked_chunk, norms_chunk,
                                        ranked_all, norms_all, m)

    for local_i, global_i in enumerate(chunk_indices):
        col_pos = int(np.searchsorted(all_indices, global_i))
        if col_pos < r_block.shape[1] and all_indices[col_pos] == global_i:
            r_block[local_i, col_pos] = np.nan
            p_block[local_i, col_pos] = np.nan

    pass_mask = (r_block >= min_thr) & (p_block <= p_thr)
    ia, ib = np.where(pass_mask)
    r_vals = r_block[ia, ib]
    p_vals = p_block[ia, ib]

    empty_cols = ["seed", "target", "r", "p", "r_group",
                  "seed_nonzero_samples", "target_nonzero_samples",
                  "nonzero_intersection", "nonzero_union"]

    if r_vals.size == 0:
        shm_data.close(); shm_ranked.close(); shm_norms.close()
        return pd.DataFrame(columns=empty_cols)

    presence_chunk = (data[chunk_indices] > 0)
    presence_all   = (data > 0)

    r_group       = np.array([max(t for t in thresholds if rv >= t) for rv in r_vals])
    seeds_out     = labels_chunk[ia]
    targets_out   = labels_all[ib]

    intersections = (presence_chunk[ia] & presence_all[ib]).sum(axis=1)
    unions        = (presence_chunk[ia] | presence_all[ib]).sum(axis=1)

    shm_data.close()
    shm_ranked.close()
    shm_norms.close()

    return pd.DataFrame({
        "seed":                   seeds_out,
        "target":                 targets_out,
        "r":                      r_vals,
        "p":                      p_vals,
        "r_group":                r_group,
        "seed_nonzero_samples":   [nonzero_map.get(s) for s in seeds_out],
        "target_nonzero_samples": [nonzero_map.get(t) for t in targets_out],
        "nonzero_intersection":   intersections,
        "nonzero_union":          unions,
    })


def main():
    parser = argparse.ArgumentParser(
        description="All-to-all Spearman correlation — filtered output only."
    )
    parser.add_argument("--input",           required=True)
    parser.add_argument("--output-filtered", required=True)
    parser.add_argument("--min-samples", type=int, default=0)
    parser.add_argument("--thresholds", default="0.6,0.7,0.8,0.9")
    parser.add_argument("--p-threshold", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--chunk-size", type=int, default=300)
    parser.add_argument("--value-suffix", default="_rpkmf",
                         help="suffix identifying per-sample value columns (default: _rpkmf)")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output_filtered)), exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    rpkmf_cols = [c for c in df.columns if c.endswith(args.value_suffix)]
    if not rpkmf_cols:
        raise SystemExit(f"No {args.value_suffix} columns found in input.")

    label_col = _label_col(df, args.value_suffix)
    for c in rpkmf_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    if args.min_samples > 0:
        presence = (df[rpkmf_cols] > 0).sum(axis=1)
        df = df[presence >= args.min_samples].reset_index(drop=True)

    labels = df[label_col].astype(str).tolist()
    data   = df[rpkmf_cols].to_numpy(dtype=float)

    nonzero_counts = (data > 0).sum(axis=1)
    nonzero_map    = dict(zip(labels, nonzero_counts.tolist()))

    n, m = data.shape
    thresholds = _parse_thresholds(args.thresholds)
    min_thr    = thresholds[0]

    empty_cols = ["seed", "target", "r", "p", "r_group",
                  "seed_nonzero_samples", "target_nonzero_samples",
                  "nonzero_intersection", "nonzero_union"]

    if n < 2 or m < 3:
        pd.DataFrame(columns=empty_cols).to_csv(args.output_filtered, sep="\t", index=False)
        return

    labels_arr = np.array(labels)
    all_idx    = np.arange(n)
    chunk      = args.chunk_size
    n_workers  = args.threads
    chunk_list = [all_idx[s:min(s + chunk, n)] for s in range(0, n, chunk)]
    n_chunks   = len(chunk_list)

    print(f"Input: {n} rows x {m} samples  (all-to-all)", file=sys.stderr, flush=True)
    print(f"chunk_size={chunk}  threads={n_workers}  total chunks={n_chunks}",
          file=sys.stderr, flush=True)

    data_c       = np.ascontiguousarray(data, dtype=np.float64)
    shm_data     = SharedMemory(create=True, size=data_c.nbytes)
    shm_data_arr = np.ndarray(data_c.shape, dtype=data_c.dtype, buffer=shm_data.buf)
    shm_data_arr[:] = data_c

    def _make_shm(arr: np.ndarray):
        arr_c = np.ascontiguousarray(arr, dtype=np.float64)
        shm   = SharedMemory(create=True, size=arr_c.nbytes)
        buf   = np.ndarray(arr_c.shape, dtype=np.float64, buffer=shm.buf)
        buf[:] = arr_c
        return shm, arr_c.shape

    print("Pre-ranking all rows ...", file=sys.stderr, flush=True)
    ranked_all, norms_all = _prerank(data)
    shm_ranked, ranked_shape = _make_shm(ranked_all)
    shm_norms,  norms_shape  = _make_shm(norms_all)

    out_path = args.output_filtered
    with open(out_path, "w") as fh:
        fh.write("\t".join(empty_cols) + "\n")

    total_written = 0
    max_inflight = n_workers * 2

    try:
        with ProcessPoolExecutor(max_workers=n_workers) as pool:
            pending = {}
            ci_next = 0
            completed = 0

            def _submit_next():
                nonlocal ci_next
                while ci_next < n_chunks and len(pending) < max_inflight:
                    chunk_indices = chunk_list[ci_next]
                    f = pool.submit(
                        _worker,
                        chunk_indices, all_idx,
                        shm_data.name, data_c.shape, data_c.dtype,
                        shm_ranked.name, ranked_shape,
                        shm_norms.name,  norms_shape,
                        m,
                        labels_arr[chunk_indices],
                        labels_arr,
                        min_thr, args.p_threshold, thresholds,
                        nonzero_map,
                    )
                    pending[f] = ci_next
                    ci_next += 1

            _submit_next()

            while pending:
                for f in as_completed(pending):
                    df_chunk = f.result()
                    del pending[f]
                    if len(df_chunk) > 0:
                        df_chunk.to_csv(out_path, sep="\t", index=False,
                                        header=False, mode="a")
                        total_written += len(df_chunk)
                    completed += 1
                    print(f"  chunk {completed}/{n_chunks} done  "
                          f"(+{len(df_chunk)} pairs, {total_written} total)",
                          file=sys.stderr, flush=True)
                    _submit_next()
                    break

    finally:
        shm_data.close(); shm_data.unlink()
        shm_ranked.close(); shm_ranked.unlink()
        shm_norms.close();  shm_norms.unlink()

    print(f"Done. {total_written} pairs passed filters -> {out_path}",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
