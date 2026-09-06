#!/usr/bin/env python3

import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing.shared_memory import SharedMemory
from typing import List

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as t_dist


# ------------------------------------------------------------------ #
# Pure functions (picklable, called in worker processes)              #
# ------------------------------------------------------------------ #

def _is_unknown(value: str) -> bool:
    v = str(value).strip()
    vl = v.lower()
    return (vl in {"unknow", "unknown", ""}
            or vl.startswith("unknown::")
            or vl.startswith("unknow::")
            or v.startswith("acc:"))


def _parse_thresholds(text: str) -> List[float]:
    vals = [float(p.strip()) for p in text.split(",") if p.strip()]
    if not vals:
        raise SystemExit("No valid thresholds provided.")
    return sorted(set(vals))


def _label_col(df: pd.DataFrame) -> str:
    for c in df.columns:
        if not c.endswith("_RPKMF"):
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
    chunk_indices: np.ndarray,   # row indices into data for this chunk (seed side)
    fixed_indices: np.ndarray,   # row indices into data for the fixed side (target)
    # shared memory descriptors for data matrix
    shm_data_name: str, data_shape: tuple, data_dtype,
    # shared memory descriptors for pre-ranked fixed side
    shm_ranked_name: str, ranked_shape: tuple,
    shm_norms_name:  str, norms_shape:  tuple,
    m: int,
    labels_chunk: np.ndarray,
    labels_fixed: np.ndarray,
    min_thr: float,
    p_thr: float,
    thresholds: List[float],
    exclude_diagonal: bool,   # True for kk/uu where seed==target array
):
    """Worker: compute one chunk, return (df_all, df_filt)."""
    # attach shared memory (read-only views)
    shm_data   = SharedMemory(name=shm_data_name,   create=False)
    shm_ranked = SharedMemory(name=shm_ranked_name, create=False)
    shm_norms  = SharedMemory(name=shm_norms_name,  create=False)

    data        = np.ndarray(data_shape,   dtype=data_dtype,  buffer=shm_data.buf)
    ranked_fixed = np.ndarray(ranked_shape, dtype=np.float64, buffer=shm_ranked.buf)
    norms_fixed  = np.ndarray(norms_shape,  dtype=np.float64, buffer=shm_norms.buf)

    ranked_chunk, norms_chunk = _prerank(data[chunk_indices])
    r_block, p_block = _spearman_block(ranked_chunk, norms_chunk,
                                        ranked_fixed, norms_fixed, m)

    if exclude_diagonal:
        for local_i, global_i in enumerate(chunk_indices):
            col_pos = int(np.searchsorted(fixed_indices, global_i))
            if col_pos < r_block.shape[1] and fixed_indices[col_pos] == global_i:
                r_block[local_i, col_pos] = np.nan
                p_block[local_i, col_pos] = np.nan

    ia, ib = np.where(np.isfinite(r_block))
    r_vals = r_block[ia, ib]
    p_vals = p_block[ia, ib]

    df_all = pd.DataFrame({
        "seed":   labels_chunk[ia],
        "target": labels_fixed[ib],
        "r":      r_vals,
        "p":      p_vals,
    })

    filt = (r_vals >= min_thr) & (p_vals <= p_thr)
    if filt.any():
        r_f = r_vals[filt]
        r_group = np.array([max(t for t in thresholds if rv >= t) for rv in r_f])
        df_filt = pd.DataFrame({
            "seed":    labels_chunk[ia[filt]],
            "target":  labels_fixed[ib[filt]],
            "r":       r_f,
            "p":       p_vals[filt],
            "r_group": r_group,
        })
    else:
        df_filt = pd.DataFrame(columns=["seed", "target", "r", "p", "r_group"])

    shm_data.close()
    shm_ranked.close()
    shm_norms.close()

    return df_all, df_filt


def _run_pair_type(
    tag: str,
    chunk_idx_list: List[np.ndarray],  # one array of row indices per chunk (seed side)
    fixed_indices: np.ndarray,          # all row indices for the fixed (target) side
    shm_data_name: str, data_shape, data_dtype,
    shm_ranked_name: str, ranked_shape,
    shm_norms_name:  str, norms_shape,
    m: int,
    labels_arr: np.ndarray,
    min_thr: float, p_thr: float, thresholds: List[float],
    exclude_diagonal: bool,
    n_workers: int,
    rows_all: list, rows_filt: list,
    n_seed_total: int,
    n_fixed_total: int,
):
    n_chunks = len(chunk_idx_list)
    print(f"{tag} pairs: {n_seed_total} × {n_fixed_total} "
          f"→ {n_chunks} chunks  ({n_workers} workers)",
          file=sys.stderr, flush=True)

    futures = {}
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        for ci, chunk_indices in enumerate(chunk_idx_list):
            f = pool.submit(
                _worker,
                chunk_indices, fixed_indices,
                shm_data_name, data_shape, data_dtype,
                shm_ranked_name, ranked_shape,
                shm_norms_name,  norms_shape,
                m,
                labels_arr[chunk_indices],
                labels_arr[fixed_indices],
                min_thr, p_thr, thresholds,
                exclude_diagonal,
            )
            futures[f] = ci

        completed = 0
        for f in as_completed(futures):
            ci = futures[f]
            df_all_chunk, df_filt_chunk = f.result()
            rows_all.append(df_all_chunk)
            rows_filt.append(df_filt_chunk)
            completed += 1
            chunk_indices = chunk_idx_list[ci]
            start = int(chunk_indices[0] - chunk_idx_list[0][0] if ci == 0
                        else chunk_idx_list[ci][0] - chunk_idx_list[0][0])
            print(f"  [{tag}] chunk {completed}/{n_chunks} done",
                  file=sys.stderr, flush=True)


def main():
    parser = argparse.ArgumentParser(
        description="Pairwise Spearman correlation on any pre-aggregated RPKMF table."
    )
    parser.add_argument("--input",           required=True)
    parser.add_argument("--output-all",      required=True)
    parser.add_argument("--output-filtered", required=True)
    parser.add_argument("--min-samples", type=int, default=0)
    parser.add_argument("--thresholds", default="0.6,0.7,0.8,0.9")
    parser.add_argument("--p-threshold", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=1,
                        help="Parallel worker processes (default: 1)")
    parser.add_argument("--pairs", default="kk,uk",
                        help="kk=known-known, uk=unknown-known, ku=known-unknown, "
                             "uu=unknown-unknown  (default: kk,uk)")
    parser.add_argument("--clr", action="store_true", default=False)
    parser.add_argument("--clr-pseudocount", type=float, default=None)
    parser.add_argument("--chunk-size", type=int, default=5000,
                        help="Rows per chunk (default: 5000)")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.output_all)),      exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.output_filtered)), exist_ok=True)

    df = pd.read_csv(args.input, sep="\t", dtype=str)
    rpkmf_cols = [c for c in df.columns if c.endswith("_RPKMF")]
    if not rpkmf_cols:
        raise SystemExit("No _RPKMF columns found in input.")

    label_col = _label_col(df)
    for c in rpkmf_cols:
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)

    if args.min_samples > 0:
        presence = (df[rpkmf_cols] > 0).sum(axis=1)
        df = df[presence >= args.min_samples].reset_index(drop=True)

    labels = df[label_col].astype(str).tolist()
    raw = df[rpkmf_cols].to_numpy(dtype=float)

    if args.clr:
        pseudo = args.clr_pseudocount
        if pseudo is None:
            nonzero = raw[raw > 0]
            pseudo = nonzero.min() / 2 if nonzero.size > 0 else 1.0
        log_shifted = np.log(raw + pseudo)
        data = log_shifted - log_shifted.mean(axis=0, keepdims=True)
    else:
        data = raw

    n, m = data.shape
    thresholds = _parse_thresholds(args.thresholds)
    min_thr = thresholds[0]

    # non-zero sample counts per virus (used to annotate output)
    nonzero_counts = (raw > 0).sum(axis=1)
    nonzero_map = dict(zip(labels, nonzero_counts.tolist()))

    empty_header_all  = ["seed", "target", "r", "p", "seed_nonzero_samples", "target_nonzero_samples"]
    empty_header_filt = ["seed", "target", "r", "p", "r_group", "seed_nonzero_samples", "target_nonzero_samples"]

    if n < 2 or m < 3:
        pd.DataFrame(columns=empty_header_all).to_csv(args.output_all, sep="\t", index=False)
        pd.DataFrame(columns=empty_header_filt).to_csv(args.output_filtered, sep="\t", index=False)
        return

    pair_types = {p.strip().lower() for p in args.pairs.split(",")}
    invalid = pair_types - {"kk", "uk", "ku", "uu"}
    if invalid:
        raise SystemExit(f"Unknown pair types: {invalid}. Choose from: kk, uk, ku, uu")

    is_unk    = np.array([_is_unknown(l) for l in labels])
    labels_arr = np.array(labels)
    n_known   = int((~is_unk).sum())
    n_unknown = int(is_unk.sum())

    print(f"Input: {n} rows ({n_known} known, {n_unknown} unknown) x {m} samples",
          file=sys.stderr, flush=True)
    print(f"Pair types: {', '.join(sorted(pair_types))}  "
          f"chunk_size={args.chunk_size}  threads={args.threads}",
          file=sys.stderr, flush=True)

    unk_idx   = np.where(is_unk)[0]
    known_idx = np.where(~is_unk)[0]
    chunk     = args.chunk_size
    n_workers = args.threads

    # ------------------------------------------------------------------ #
    # Put data matrix in shared memory so workers don't copy it           #
    # ------------------------------------------------------------------ #
    data_c = np.ascontiguousarray(data, dtype=np.float64)
    shm_data = SharedMemory(create=True, size=data_c.nbytes)
    shm_data_arr = np.ndarray(data_c.shape, dtype=data_c.dtype, buffer=shm_data.buf)
    shm_data_arr[:] = data_c

    def _make_shm(arr: np.ndarray):
        arr_c = np.ascontiguousarray(arr, dtype=np.float64)
        shm = SharedMemory(create=True, size=arr_c.nbytes)
        buf = np.ndarray(arr_c.shape, dtype=np.float64, buffer=shm.buf)
        buf[:] = arr_c
        return shm, arr_c.shape

    need_kk = "kk" in pair_types
    need_uk = "uk" in pair_types
    need_ku = "ku" in pair_types
    need_uu = "uu" in pair_types

    # pre-rank fixed sides and put in shared memory
    shm_ranked_known = shm_norms_known = None
    shm_ranked_unk   = shm_norms_unk   = None
    ranked_known_shape = norms_known_shape = None
    ranked_unk_shape   = norms_unk_shape   = None

    if (need_kk or need_uk or need_ku) and len(known_idx) > 0:
        print("Pre-ranking known side ...", file=sys.stderr, flush=True)
        ranked_known, norms_known = _prerank(data[known_idx])
        shm_ranked_known, ranked_known_shape = _make_shm(ranked_known)
        shm_norms_known,  norms_known_shape  = _make_shm(norms_known)

    if (need_uu or need_ku or need_uk) and len(unk_idx) > 0:
        print("Pre-ranking unknown side ...", file=sys.stderr, flush=True)
        ranked_unk, norms_unk = _prerank(data[unk_idx])
        shm_ranked_unk, ranked_unk_shape = _make_shm(ranked_unk)
        shm_norms_unk,  norms_unk_shape  = _make_shm(norms_unk)

    rows_all  = []
    rows_filt = []

    def _chunks_of(idx):
        return [idx[s:min(s+chunk, len(idx))]
                for s in range(0, len(idx), chunk)]

    try:
        if need_kk and len(known_idx) > 0:
            _run_pair_type(
                "kk", _chunks_of(known_idx), known_idx,
                shm_data.name, data_c.shape, data_c.dtype,
                shm_ranked_known.name, ranked_known_shape,
                shm_norms_known.name,  norms_known_shape,
                m, labels_arr, min_thr, args.p_threshold, thresholds,
                exclude_diagonal=True,
                n_workers=n_workers,
                rows_all=rows_all, rows_filt=rows_filt,
                n_seed_total=len(known_idx), n_fixed_total=len(known_idx),
            )

        if need_uu and len(unk_idx) > 0:
            _run_pair_type(
                "uu", _chunks_of(unk_idx), unk_idx,
                shm_data.name, data_c.shape, data_c.dtype,
                shm_ranked_unk.name, ranked_unk_shape,
                shm_norms_unk.name,  norms_unk_shape,
                m, labels_arr, min_thr, args.p_threshold, thresholds,
                exclude_diagonal=True,
                n_workers=n_workers,
                rows_all=rows_all, rows_filt=rows_filt,
                n_seed_total=len(unk_idx), n_fixed_total=len(unk_idx),
            )

        if need_uk and len(unk_idx) > 0 and len(known_idx) > 0:
            _run_pair_type(
                "uk", _chunks_of(unk_idx), known_idx,
                shm_data.name, data_c.shape, data_c.dtype,
                shm_ranked_known.name, ranked_known_shape,
                shm_norms_known.name,  norms_known_shape,
                m, labels_arr, min_thr, args.p_threshold, thresholds,
                exclude_diagonal=False,
                n_workers=n_workers,
                rows_all=rows_all, rows_filt=rows_filt,
                n_seed_total=len(unk_idx), n_fixed_total=len(known_idx),
            )

        if need_ku and len(known_idx) > 0 and len(unk_idx) > 0:
            _run_pair_type(
                "ku", _chunks_of(known_idx), unk_idx,
                shm_data.name, data_c.shape, data_c.dtype,
                shm_ranked_unk.name, ranked_unk_shape,
                shm_norms_unk.name,  norms_unk_shape,
                m, labels_arr, min_thr, args.p_threshold, thresholds,
                exclude_diagonal=False,
                n_workers=n_workers,
                rows_all=rows_all, rows_filt=rows_filt,
                n_seed_total=len(known_idx), n_fixed_total=len(unk_idx),
            )

    finally:
        # always release shared memory
        shm_data.close(); shm_data.unlink()
        for shm in [shm_ranked_known, shm_norms_known,
                    shm_ranked_unk,   shm_norms_unk]:
            if shm is not None:
                shm.close(); shm.unlink()

    print("Writing output ...", file=sys.stderr, flush=True)

    def _add_nonzero_cols(df):
        df["seed_nonzero_samples"]   = df["seed"].map(nonzero_map)
        df["target_nonzero_samples"] = df["target"].map(nonzero_map)
        return df

    df_all = pd.concat(rows_all, ignore_index=True) if rows_all else pd.DataFrame(columns=empty_header_all)
    df_all = _add_nonzero_cols(df_all)
    df_all.to_csv(args.output_all, sep="\t", index=False)

    df_filt = pd.concat(rows_filt, ignore_index=True) if rows_filt else pd.DataFrame(columns=empty_header_filt)
    df_filt = _add_nonzero_cols(df_filt)
    df_filt.to_csv(args.output_filtered, sep="\t", index=False)

    print(f"Done. {len(df_all)} pairs total, {len(df_filt)} passed filters.",
          file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
