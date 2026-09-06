# Line-by-Line Code Explanation: esviritu_subspecies_correlation_v2.py

This document explains every single line of the script in plain language.
No coding background is assumed. Each block of code is shown first, then
explained line by line.

---

## Block 1 — The shebang line (line 1)

```python
#!/usr/bin/env python3
```

**Line 1:** This is called a "shebang". It tells the operating system: "run
this file using Python 3". It is only needed when you execute the script
directly from the terminal (e.g. `./script.py`). When called via
`python script.py`, this line is ignored. It has no effect on the logic.

---

## Block 2 — Imports (lines 3–12)

```python
import argparse
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing.shared_memory import SharedMemory
from typing import List

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as t_dist
```

**Line 3 — `import argparse`:**
`argparse` is a built-in Python library that reads command-line arguments.
When Snakemake runs the script with flags like `--input file.tsv --threads 16`,
`argparse` parses those flags and makes their values available inside the script.

**Line 4 — `import os`:**
`os` gives access to operating system functions. Here it is used only to
create output directories (`os.makedirs`).

**Line 5 — `import sys`:**
`sys` gives access to system-level features. Here it is used only for
`sys.stderr` — the "error/log" output stream. Progress messages are printed
there so they appear in the `.err` log file and not in the output data file.

**Line 6 — `from concurrent.futures import ProcessPoolExecutor, as_completed`:**
This imports the tools for running code in parallel across multiple CPU cores.
- `ProcessPoolExecutor` — a "pool" of worker processes. You tell it how many
  workers (CPUs) to use, and it manages distributing work to them.
- `as_completed` — an iterator that yields results as each worker finishes,
  rather than waiting for all workers to finish before processing any result.

**Line 7 — `from multiprocessing.shared_memory import SharedMemory`:**
`SharedMemory` creates a block of RAM that multiple processes can all read
simultaneously without copying data between them. This is critical here because
the main data matrix (~62 MB) would otherwise be copied to every worker process
— wasting time and memory.

**Line 8 — `from typing import List`:**
`List` is used in type hints like `List[float]`, which just means "a list
containing float numbers". This is documentation only — Python does not enforce
it at runtime. It helps human readers understand what a function expects.

**Line 10 — `import numpy as np`:**
NumPy is the core scientific computing library. It provides fast multi-dimensional
arrays and mathematical operations (matrix multiply, ranking, norms). The alias
`np` means you write `np.array(...)` instead of `numpy.array(...)`.

**Line 11 — `import pandas as pd`:**
Pandas provides the `DataFrame` — a table with named columns, like a spreadsheet.
Used here to read the input TSV file and write the output TSV files. The alias
`pd` is standard convention.

**Line 12 — `from scipy.stats import rankdata, t as t_dist`:**
- `rankdata` — converts a list of numbers into their ranks (1st, 2nd, 3rd...).
  For example `rankdata([10, 30, 20])` → `[1.0, 3.0, 2.0]`.
- `t as t_dist` — imports the t-distribution under the name `t_dist`. Used to
  convert t-statistics into p-values. The rename avoids confusion with the
  variable name `t` used inside `_spearman_block`.

---

## Block 3 — `_is_unknown()` function (lines 19–25)

```python
def _is_unknown(value: str) -> bool:
    v = str(value).strip()
    vl = v.lower()
    return (vl in {"unknow", "unknown", ""}
            or vl.startswith("unknown::")
            or vl.startswith("unknow::")
            or v.startswith("acc:"))
```

**What this function does:** Given a virus label (its name/ID), decides whether
it represents an "unknown" virus (returns `True`) or a "known" virus (returns
`False`).

**Line `def _is_unknown(value: str) -> bool:`:**
Defines the function. It takes one input (`value`, a string) and returns a
boolean (`True` or `False`). The leading underscore in `_is_unknown` is a
Python convention meaning "this is an internal helper, not meant to be called
from outside this script".

**Line `v = str(value).strip()`:**
- `str(value)` — converts the input to a string, just in case it came in as
  something else (e.g. a number or NaN from pandas).
- `.strip()` — removes any leading or trailing whitespace. For example
  `"  unknown  "` becomes `"unknown"`. Without this, a label with an accidental
  space would not be detected.

**Line `vl = v.lower()`:**
Creates a lowercase copy of the label. For example `"UNKNOWN"` → `"unknown"`.
This makes the check case-insensitive, so `"Unknown"`, `"UNKNOWN"`, `"unknown"`
all behave the same.

**Line `return (vl in {"unknow", "unknown", ""})`:**
Returns `True` if the lowercase label is exactly one of these three values:
- `"unknown"` — the standard unknown label
- `"unknow"` — a common typo found in some EsViritu outputs
- `""` — empty string, meaning no label was assigned

`{"unknow", "unknown", ""}` is a Python **set** (curly braces). Checking
membership in a set (`vl in {...}`) is faster than checking in a list.

**Line `or vl.startswith("unknown::")`:**
Also returns `True` if the label starts with `"unknown::"`. These labels are
created in the merge step when a virus could not be assigned to a named
subspecies — it gets the label `"unknown::acc:SampleName_contigID"` to make
it uniquely identifiable while flagging it as unknown.

**Line `or vl.startswith("unknow::")`:**
Same as above but for the typo variant `"unknow::"`.

**Line `or v.startswith("acc:")`:**
Note: this uses `v` (original case), not `vl` (lowercase). Labels starting
with `"acc:"` are custom-assembled contigs from your own sequencing data — they
are not in any reference database and so are treated as unknown at assembly
level. Example: `"acc:WanZhou_N_5_0000001133"`.

---

## Block 4 — `_parse_thresholds()` function (lines 28–32)

```python
def _parse_thresholds(text: str) -> List[float]:
    vals = [float(p.strip()) for p in text.split(",") if p.strip()]
    if not vals:
        raise SystemExit("No valid thresholds provided.")
    return sorted(set(vals))
```

**What this function does:** Converts the string `"0.6,0.7,0.8,0.9"` (as
passed on the command line) into the Python list `[0.6, 0.7, 0.8, 0.9]`.

**Line `vals = [float(p.strip()) for p in text.split(",") if p.strip()]`:**
This is a "list comprehension" — a compact way to build a list.
- `text.split(",")` — splits the string at every comma: `["0.6", "0.7", "0.8", "0.9"]`
- `for p in ...` — loops over each piece
- `p.strip()` — removes whitespace from each piece
- `if p.strip()` — skips empty pieces (e.g. a trailing comma would produce `""`)
- `float(p.strip())` — converts the string `"0.6"` to the number `0.6`

**Line `if not vals:`:**
If the list is empty (e.g. the user passed `--thresholds ""`), raise an error
and exit immediately with a helpful message.

**Line `return sorted(set(vals))`:**
- `set(vals)` — removes any duplicate values (e.g. `"0.6,0.6,0.8"` → `{0.6, 0.8}`)
- `sorted(...)` — puts them in ascending order

---

## Block 5 — `_label_col()` function (lines 35–39)

```python
def _label_col(df: pd.DataFrame) -> str:
    for c in df.columns:
        if not c.endswith("_RPKMF"):
            return c
    raise SystemExit("No label column found in input.")
```

**What this function does:** Finds the virus name column in the input table.

The input table looks like this:

```
Assembly          Sample1_RPKMF  Sample2_RPKMF  Sample3_RPKMF
GCF_000849225.1   12.4           0.0            5.1
acc:WanZhou_N_5   0.0            3.2            0.0
```

The first column (`Assembly`) is the label — the virus name. All other columns
end with `_RPKMF` and contain numbers.

**Line `for c in df.columns:`:**
Loops through each column name in the table, in order from left to right.

**Line `if not c.endswith("_RPKMF"): return c`:**
As soon as it finds a column that does NOT end with `"_RPKMF"`, it returns that
column name. Because the label column is always first, this returns `"Assembly"`
(or `"subspecies"` or `"species"` depending on which input file is used).

**Line `raise SystemExit("No label column found in input.")`:**
If somehow every column ends with `_RPKMF` (which should never happen), exit
with an error message.

---

## Block 6 — `_prerank()` function (lines 42–46)

```python
def _prerank(data: np.ndarray):
    ranked = np.apply_along_axis(rankdata, 1, data).astype(float)
    ranked -= ranked.mean(axis=1, keepdims=True)
    norms = np.linalg.norm(ranked, axis=1)
    return ranked, norms
```

**What this function does:** Transforms the raw RPKMF matrix into centred rank
vectors, which are needed to compute Spearman correlation via matrix
multiplication.

**Why ranks?** Spearman correlation measures whether two variables move
*together in the same direction*, regardless of the actual numbers. By
converting values to ranks (1st, 2nd, 3rd...), we focus only on the ordering
pattern.

**Line `ranked = np.apply_along_axis(rankdata, 1, data).astype(float)`:**
- `data` is a 2D matrix: rows = viruses, columns = samples
- `np.apply_along_axis(rankdata, 1, data)` applies the `rankdata` function to
  each row (axis=1 means "along each row"). For each virus, it converts its
  RPKMF values across all samples into ranks.
  - Example: virus with values `[0, 0, 5, 2, 8]` across 5 samples gets ranks
    `[1.5, 1.5, 3, 2, ... ]` (ties at 0 get the average rank 1.5)
- `.astype(float)` converts the result to decimal numbers (needed for the next step)

**Line `ranked -= ranked.mean(axis=1, keepdims=True)`:**
- `ranked.mean(axis=1, keepdims=True)` computes the mean rank for each virus
  (mean across its row). `keepdims=True` keeps the result as a column vector
  so it can be subtracted properly.
- `ranked -=` subtracts the mean from every value in the row, "centring" each
  row around zero.
- Why centre? The Spearman formula (as a dot product) requires centred vectors
  to give correct results between -1 and +1.
- Example: ranks `[1.5, 1.5, 3, 4, 5]`, mean = 3 → centred: `[-1.5, -1.5, 0, 1, 2]`

**Line `norms = np.linalg.norm(ranked, axis=1)`:**
- Computes the "length" (Euclidean norm) of each centred rank vector.
- For a vector `[a, b, c, ...]`, the norm = `sqrt(a² + b² + c² + ...)`
- This is used as the denominator when computing correlation (to normalise
  the result to be between -1 and +1).

**Line `return ranked, norms`:**
Returns both the centred rank matrix (shape: n_viruses × n_samples) and the
norms (shape: n_viruses).

---

## Block 7 — `_spearman_block()` function (lines 49–55)

```python
def _spearman_block(ranked_a, norms_a, ranked_b, norms_b, m):
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (ranked_a @ ranked_b.T) / np.outer(norms_a, norms_b)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = r * np.sqrt((m - 2) / (1 - r ** 2))
    p = 2 * t_dist.sf(np.abs(t), df=m - 2)
    return r, p
```

**What this function does:** Given two groups of pre-ranked viruses (A and B),
computes the Spearman correlation and p-value for every possible (A_i, B_k) pair
simultaneously using matrix operations.

**Line `with np.errstate(invalid="ignore", divide="ignore"):`:**
Suppresses NumPy's warnings about mathematical edge cases (like dividing by
zero or taking the square root of a negative number). These can happen for
viruses that are all-zero or perfectly identical, but the result (NaN or Inf)
is handled downstream. Without this, the terminal would be flooded with warnings.

**Line `r = (ranked_a @ ranked_b.T) / np.outer(norms_a, norms_b)`:**
This computes the entire correlation matrix in one operation.
- `ranked_a @ ranked_b.T` — matrix multiplication of group A's rank vectors
  by group B's rank vectors (transposed). The result is a matrix where entry
  [i, k] is the dot product of virus A_i's ranks with virus B_k's ranks.
  This is equivalent to computing the numerator of the Pearson correlation
  for each pair.
- `np.outer(norms_a, norms_b)` — creates a matrix where entry [i, k] is
  `norm_A_i × norm_B_k`. This is the denominator.
- Dividing gives the Spearman correlation r for every pair.
- Result shape: (n_A × n_B)

**Line `t = r * np.sqrt((m - 2) / (1 - r ** 2))`:**
Converts each correlation r to a t-statistic.
- `m` is the number of samples
- `r ** 2` is r squared
- `1 - r**2` represents "unexplained variance"
- The formula is the standard t-test for a correlation coefficient
- Applied to every cell of the matrix simultaneously

**Line `p = 2 * t_dist.sf(np.abs(t), df=m - 2)`:**
Converts t-statistics to p-values.
- `np.abs(t)` — takes absolute value (we don't care about direction for the
  p-value)
- `t_dist.sf(x, df=m-2)` — the "survival function" of the t-distribution:
  probability of observing a t-value more extreme than x. This gives the
  one-tailed p-value.
- `2 *` — makes it two-tailed (we care about both positive and negative
  correlations)
- Result: a matrix of p-values, same shape as r

**Line `return r, p`:**
Returns both matrices.

---

## Block 8 — `_worker()` function (lines 58–124)

This is the function that each parallel CPU worker runs. It receives one chunk
of virus rows and processes them against the fixed target side.

### Setup lines (76–86)

```python
shm_data   = SharedMemory(name=shm_data_name,   create=False)
shm_ranked = SharedMemory(name=shm_ranked_name, create=False)
shm_norms  = SharedMemory(name=shm_norms_name,  create=False)

data         = np.ndarray(data_shape,   dtype=data_dtype,  buffer=shm_data.buf)
ranked_fixed = np.ndarray(ranked_shape, dtype=np.float64, buffer=shm_ranked.buf)
norms_fixed  = np.ndarray(norms_shape,  dtype=np.float64, buffer=shm_norms.buf)
```

**Lines `shm_data = SharedMemory(name=..., create=False)`:**
Each worker process attaches to the shared memory blocks that the main process
already created. `create=False` means "attach to existing block, don't create
a new one". The `name` is a unique identifier (like a filename) for the block.

**Lines `data = np.ndarray(..., buffer=shm_data.buf)`:**
Creates a NumPy array that points directly into the shared memory block — no
copying. When the worker reads `data[i]`, it reads from the shared block.
`data_shape` and `data_dtype` tell NumPy how to interpret the raw bytes (e.g.
"this is a 161760 × 48 matrix of 64-bit floats").

### Core computation (84–86)

```python
ranked_chunk, norms_chunk = _prerank(data[chunk_indices])
r_block, p_block = _spearman_block(ranked_chunk, norms_chunk,
                                    ranked_fixed, norms_fixed, m)
```

**Line `ranked_chunk, norms_chunk = _prerank(data[chunk_indices])`:**
- `data[chunk_indices]` selects only the rows for this chunk from the full
  data matrix. For example, `chunk_indices = [5000, 5001, ..., 9999]` selects
  rows 5000 to 9999.
- Calls `_prerank()` to rank and centre those rows.

**Lines `r_block, p_block = _spearman_block(...)`:**
Computes the correlation and p-value between every virus in this chunk (seed)
and every virus in the fixed target group. Returns two matrices of shape
(chunk_size × n_fixed).

### Diagonal exclusion (88–93)

```python
if exclude_diagonal:
    for local_i, global_i in enumerate(chunk_indices):
        col_pos = int(np.searchsorted(fixed_indices, global_i))
        if col_pos < r_block.shape[1] and fixed_indices[col_pos] == global_i:
            r_block[local_i, col_pos] = np.nan
            p_block[local_i, col_pos] = np.nan
```

**Line `if exclude_diagonal:`:**
Only applies for `kk` and `uu` pair types, where the seed and target groups are
the same. In that case, virus i would be correlated with itself — giving r = 1.0
always, which is not informative and should be excluded.

**Line `for local_i, global_i in enumerate(chunk_indices):`:**
Loops over each virus in this chunk. `local_i` is the row position within the
chunk (0, 1, 2...), `global_i` is the row index in the full data matrix.

**Line `col_pos = int(np.searchsorted(fixed_indices, global_i))`:**
Finds where `global_i` sits in the sorted `fixed_indices` array. This gives
the column position in `r_block` that corresponds to the same virus.

**Lines `r_block[local_i, col_pos] = np.nan`:**
Sets the self-correlation cell to NaN (Not a Number). NaN values are
automatically excluded when the results are flattened in the next step.

### Flattening results (95–104)

```python
ia, ib = np.where(np.isfinite(r_block))
r_vals = r_block[ia, ib]
p_vals = p_block[ia, ib]

df_all = pd.DataFrame({
    "seed":   labels_chunk[ia],
    "target": labels_fixed[ib],
    "r":      r_vals,
    "p":      p_vals,
})
```

**Line `ia, ib = np.where(np.isfinite(r_block))`:**
`np.isfinite` returns True for every cell that is a valid number (not NaN, not
Inf). `np.where` gives the row indices (`ia`) and column indices (`ib`) of all
True positions. This effectively finds all valid virus pairs in the block.

**Lines `r_vals = r_block[ia, ib]`:**
Uses the indices to extract the r values (and p values) for all valid pairs
into flat 1D arrays.

**Lines `df_all = pd.DataFrame({...})`:**
Creates a table with columns `seed`, `target`, `r`, `p` — one row per pair.
- `labels_chunk[ia]` — gets the seed virus name for each pair using its row index
- `labels_fixed[ib]` — gets the target virus name for each pair using its column index

### Filtering (106–118)

```python
filt = (r_vals >= min_thr) & (p_vals <= p_thr)
if filt.any():
    r_f = r_vals[filt]
    r_group = np.array([max(t for t in thresholds if rv >= t) for rv in r_f])
    df_filt = pd.DataFrame({...})
else:
    df_filt = pd.DataFrame(columns=["seed", "target", "r", "p", "r_group"])
```

**Line `filt = (r_vals >= min_thr) & (p_vals <= p_thr)`:**
Creates a boolean array — `True` for pairs that pass both thresholds:
- `r_vals >= min_thr` — correlation is strong enough (e.g. r ≥ 0.6)
- `p_vals <= p_thr` — result is statistically significant (e.g. p ≤ 0.05)
- `&` — both conditions must be true simultaneously

**Line `if filt.any():`:**
Checks if at least one pair passed the filter. If none passed, skip building
the filtered DataFrame (saves time).

**Line `r_group = np.array([max(t for t in thresholds if rv >= t) for rv in r_f])`:**
For each filtered pair, finds the highest threshold it qualifies for.
For example, if r = 0.85 and thresholds = [0.6, 0.7, 0.8, 0.9]:
- 0.85 ≥ 0.6 ✓, 0.85 ≥ 0.7 ✓, 0.85 ≥ 0.8 ✓, 0.85 ≥ 0.9 ✗
- So `max(...)` = 0.8 → this pair goes in the "0.8" group

**Lines `else: df_filt = pd.DataFrame(columns=[...])`:**
If no pairs passed the filter, returns an empty table with the correct columns
(so the output file still has proper headers).

### Cleanup (120–122)

```python
shm_data.close()
shm_ranked.close()
shm_norms.close()
```

**Each `.close()` line:**
Detaches this worker process from the shared memory blocks. Does not delete the
blocks — just says "I'm done using them". The main process is responsible for
deletion (`unlink()`).

---

## Block 9 — `_run_pair_type()` function (lines 127–176)

```python
def _run_pair_type(tag, chunk_idx_list, fixed_indices, ...):
    n_chunks = len(chunk_idx_list)
    print(f"{tag} pairs: {n_seed_total} × {n_fixed_total} "
          f"→ {n_chunks} chunks  ({n_workers} workers)",
          file=sys.stderr, flush=True)

    futures = {}
    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        for ci, chunk_indices in enumerate(chunk_idx_list):
            f = pool.submit(_worker, chunk_indices, fixed_indices, ...)
            futures[f] = ci

        completed = 0
        for f in as_completed(futures):
            ci = futures[f]
            df_all_chunk, df_filt_chunk = f.result()
            rows_all.append(df_all_chunk)
            rows_filt.append(df_filt_chunk)
            completed += 1
            print(f"  [{tag}] chunk {completed}/{n_chunks} done",
                  file=sys.stderr, flush=True)
```

**Line `n_chunks = len(chunk_idx_list)`:**
How many chunks the seed side was split into.

**Lines `print(f"{tag} pairs: ...")`:**
Prints a summary line to stderr at the start, e.g.:
`uk pairs: 161387 × 373 → 33 chunks  (16 workers)`

**Line `futures = {}`:**
A dictionary that will map each submitted task (called a "future") to its chunk
number `ci`. Needed later to report which chunk finished.

**Line `with ProcessPoolExecutor(max_workers=n_workers) as pool:`:**
Creates the worker pool with `n_workers` parallel processes (e.g. 16). The
`with` statement ensures the pool is properly shut down when done, even if an
error occurs.

**Lines `for ci, chunk_indices in enumerate(chunk_idx_list):`:**
Loops over each chunk. `ci` is the chunk number (0, 1, 2...), `chunk_indices`
is the array of row indices for that chunk.

**Line `f = pool.submit(_worker, chunk_indices, fixed_indices, ...)`:**
Submits `_worker` to run in the pool with these arguments. This does NOT wait
for the worker to finish — it immediately returns a "future" object `f`, which
is a promise that the result will be available later. All chunks are submitted
immediately, and the pool decides when to actually run them based on available
workers.

**Line `futures[f] = ci`:**
Records which chunk number belongs to this future, so when it finishes we
know which chunk it was.

**Line `for f in as_completed(futures):`:**
Waits for results to come in, processing each one as soon as it finishes
(in any order). This is more efficient than waiting for all workers to finish
before processing any result.

**Line `df_all_chunk, df_filt_chunk = f.result()`:**
Gets the return value from the finished worker — two DataFrames.

**Lines `rows_all.append(df_all_chunk)`:**
Adds the chunk's results to the running list. All chunks' results will be
concatenated at the end.

**Lines `completed += 1` and `print(...)`:**
Counts completed chunks and prints progress to stderr. With 33 chunks you will
see 33 lines like `[uk] chunk 5/33 done`.

---

## Block 10 — `main()` function (lines 179–377)

This is where the script starts executing. Everything above is function
definitions — they do nothing until called from `main()`.

### Argument parsing (lines 180–198)

```python
parser = argparse.ArgumentParser(...)
parser.add_argument("--input",           required=True)
parser.add_argument("--output-all",      required=True)
parser.add_argument("--output-filtered", required=True)
parser.add_argument("--min-samples", type=int, default=0)
parser.add_argument("--thresholds", default="0.6,0.7,0.8,0.9")
parser.add_argument("--p-threshold", type=float, default=0.05)
parser.add_argument("--threads", type=int, default=1)
parser.add_argument("--pairs", default="kk,uk")
parser.add_argument("--clr", action="store_true", default=False)
parser.add_argument("--clr-pseudocount", type=float, default=None)
parser.add_argument("--chunk-size", type=int, default=5000)
args = parser.parse_args()
```

**Line `parser = argparse.ArgumentParser(...)`:**
Creates a new argument parser.

**Each `parser.add_argument(...)` line:**
Registers one command-line flag. Key parts:
- `required=True` — Snakemake must provide this; the script will crash without it
- `type=int` / `type=float` — automatically converts the string from the
  command line to a number
- `default=...` — value used if the flag is not provided
- `action="store_true"` — for `--clr`, no value is needed; just writing
  `--clr` sets it to `True`

**Line `args = parser.parse_args()`:**
Actually reads the command-line arguments and stores them. After this line,
`args.input` contains the input file path, `args.threads` contains the number
of threads, etc.

### Creating output directories (lines 200–201)

```python
os.makedirs(os.path.dirname(os.path.abspath(args.output_all)),      exist_ok=True)
os.makedirs(os.path.dirname(os.path.abspath(args.output_filtered)), exist_ok=True)
```

**Line `os.path.abspath(args.output_all)`:**
Converts the output path to an absolute path (e.g. relative `result/05/.../file.tsv`
becomes `/scratch/user/project/result/05/.../file.tsv`). This ensures the
directory is created in the right place regardless of the current working
directory.

**Line `os.path.dirname(...)`:**
Extracts the directory part of the path, removing the filename.
e.g. `/scratch/.../1_spearman/assembly.spearman.all.tsv` → `/scratch/.../1_spearman`

**Line `os.makedirs(..., exist_ok=True)`:**
Creates the directory and any missing parent directories. `exist_ok=True` means
"if the directory already exists, don't raise an error — just continue".
This fixed the `OSError: Cannot save file into a non-existent directory` error
seen in an earlier run.

### Loading the input table (lines 203–214)

```python
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
```

**Line `df = pd.read_csv(args.input, sep="\t", dtype=str)`:**
Reads the TSV file into a DataFrame. `sep="\t"` means tab-separated (not
comma). `dtype=str` means read every column as text initially — this prevents
pandas from misinterpreting values like `"0001"` as the number `1`.

**Line `rpkmf_cols = [c for c in df.columns if c.endswith("_RPKMF")]`:**
Builds a list of column names that end with `"_RPKMF"`. These are the sample
data columns (e.g. `GouBa_N_1_RPKMF`, `WanZhou_N_5_RPKMF`...).

**Line `if not rpkmf_cols: raise SystemExit(...)`:**
If no such columns exist, the input file is wrong — exit immediately.

**Line `label_col = _label_col(df)`:**
Calls the function defined earlier to find the name column (first non-RPKMF column).

**Lines `df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0)`:**
Converts each RPKMF column from text to numbers.
- `pd.to_numeric(..., errors="coerce")` — converts text to number; any value
  that cannot be converted (e.g. `"N/A"`) becomes NaN
- `.fillna(0)` — replaces NaN with 0

**Lines `presence = (df[rpkmf_cols] > 0).sum(axis=1)`:**
For each virus (row), counts how many samples have a non-zero RPKMF.
- `df[rpkmf_cols] > 0` — True/False table: True where RPKMF > 0
- `.sum(axis=1)` — sums each row (True=1, False=0), giving count of positive samples

**Line `df = df[presence >= args.min_samples].reset_index(drop=True)`:**
Keeps only viruses detected in at least `min_samples` samples. `reset_index`
renumbers the rows from 0 after removing some.

### Normalization (lines 219–227)

```python
if args.clr:
    pseudo = args.clr_pseudocount
    if pseudo is None:
        nonzero = raw[raw > 0]
        pseudo = nonzero.min() / 2 if nonzero.size > 0 else 1.0
    log_shifted = np.log(raw + pseudo)
    data = log_shifted - log_shifted.mean(axis=0, keepdims=True)
else:
    data = raw
```

**Line `if args.clr:`:**
Only runs the CLR block if `--clr` was passed on the command line.

**Lines `pseudo = args.clr_pseudocount` / `if pseudo is None:`:**
Uses the user-specified pseudocount, or automatically computes one.

**Line `nonzero = raw[raw > 0]`:**
Selects all non-zero values from the entire matrix.

**Line `pseudo = nonzero.min() / 2`:**
Sets pseudocount to half the smallest non-zero value. This is the smallest
possible pseudocount that has a minimal effect on the data while making all
values strictly positive (needed for log).

**Line `log_shifted = np.log(raw + pseudo)`:**
Takes the natural logarithm of each value after adding the pseudocount.
`np.log(0)` would give -Inf, which is why we add `pseudo` first.

**Line `data = log_shifted - log_shifted.mean(axis=0, keepdims=True)`:**
- `log_shifted.mean(axis=0, keepdims=True)` — computes the mean of each column
  (each sample), keeping the result as a row vector so it broadcasts correctly
- Subtracts the per-sample mean from every value in that sample's column
- This centres each sample's log-values around zero — the CLR transformation

**Line `else: data = raw`:**
Without CLR, use raw RPKMF values directly. Spearman correlation uses only
ranks, so any monotone transformation (like log) would give identical results.

### Safety check (lines 236–239)

```python
if n < 2 or m < 3:
    pd.DataFrame(columns=empty_header_all).to_csv(...)
    pd.DataFrame(columns=empty_header_filt).to_csv(...)
    return
```

**Line `if n < 2 or m < 3:`:**
Spearman correlation requires at least 2 viruses (you need pairs) and at least
3 samples (with 2 samples, every non-zero correlation would be ±1, which is
meaningless). If the data is too small, write empty output files and stop.

### Classifying known vs unknown (lines 246–249)

```python
is_unk    = np.array([_is_unknown(l) for l in labels])
labels_arr = np.array(labels)
n_known   = int((~is_unk).sum())
n_unknown = int(is_unk.sum())
```

**Line `is_unk = np.array([_is_unknown(l) for l in labels])`:**
Calls `_is_unknown()` for every virus label and stores the results as a NumPy
boolean array. E.g. `[False, False, True, True, True, ...]` where `True` = unknown.

**Line `labels_arr = np.array(labels)`:**
Converts the list of labels to a NumPy array so it can be indexed by integer
arrays (needed later when extracting seed/target names).

**Lines `n_known = int((~is_unk).sum())`:**
- `~is_unk` inverts the boolean array (True → False, False → True)
- `.sum()` counts the True values = number of known viruses
- `int(...)` converts from NumPy integer to plain Python integer (for clean printing)

### Setting up shared memory (lines 265–298)

```python
data_c = np.ascontiguousarray(data, dtype=np.float64)
shm_data = SharedMemory(create=True, size=data_c.nbytes)
shm_data_arr = np.ndarray(data_c.shape, dtype=data_c.dtype, buffer=shm_data.buf)
shm_data_arr[:] = data_c
```

**Line `data_c = np.ascontiguousarray(data, dtype=np.float64)`:**
Makes a contiguous copy of the data matrix in float64 format.
"Contiguous" means the data is stored in a single unbroken block of memory
(required for shared memory). Float64 = 64-bit decimal numbers (8 bytes each).

**Line `shm_data = SharedMemory(create=True, size=data_c.nbytes)`:**
Creates a shared memory block. `data_c.nbytes` is the total size in bytes
(e.g. 161760 × 48 × 8 = ~62 MB). The OS allocates this as a named block that
other processes can attach to.

**Line `shm_data_arr = np.ndarray(..., buffer=shm_data.buf)`:**
Creates a NumPy array that views (not copies) the shared memory block.

**Line `shm_data_arr[:] = data_c`:**
Copies the actual data into shared memory. After this, `data_c` in the main
process and `shm_data_arr` point to the same bytes.

### Running pair types (lines 307–358)

```python
def _chunks_of(idx):
    return [idx[s:min(s+chunk, len(idx))]
            for s in range(0, len(idx), chunk)]
```

**Line `def _chunks_of(idx):`:**
A small helper function defined inside `main()`. Takes an array of indices and
splits it into sub-arrays of size `chunk`.

**Line `return [idx[s:min(s+chunk, len(idx))] for s in range(0, len(idx), chunk)]`:**
- `range(0, len(idx), chunk)` — generates start positions: 0, 5000, 10000, ...
- `idx[s:min(s+chunk, len(idx))]` — slices the array from `s` to `s+chunk`
  (or end of array if the last chunk is smaller)
- Returns a list of arrays, one per chunk

**Lines calling `_run_pair_type(...)` for each pair type:**
Each call computes all pairs for one type (kk, uu, uk, ku), using the
corresponding seed and fixed index arrays and shared memory handles.
`exclude_diagonal=True` is passed for kk and uu (same group on both sides),
`False` for uk and ku (different groups).

### Shared memory cleanup (lines 360–366)

```python
finally:
    shm_data.close(); shm_data.unlink()
    for shm in [shm_ranked_known, shm_norms_known,
                shm_ranked_unk,   shm_norms_unk]:
        if shm is not None:
            shm.close(); shm.unlink()
```

**`finally:` block:**
Code inside `finally` runs **always** — whether the script succeeds, crashes,
or is killed. This guarantees shared memory is always cleaned up.

**`.close()`:**
Detaches the main process from the shared memory block. The block still exists.

**`.unlink()`:**
Deletes the shared memory block from the operating system. Without this, the
memory block would persist after the script ends — consuming RAM until the
machine reboots or the block is manually deleted.

**`if shm is not None:`:**
Some blocks (e.g. `shm_ranked_unk`) are only created if the `uu` or `ku` pair
types were requested. This check avoids trying to unlink a block that was never created.

### Writing output (lines 368–374)

```python
df_all = pd.concat(rows_all, ignore_index=True) if rows_all else pd.DataFrame(columns=empty_header_all)
df_all.to_csv(args.output_all, sep="\t", index=False)

df_filt = pd.concat(rows_filt, ignore_index=True) if rows_filt else pd.DataFrame(columns=empty_header_filt)
df_filt.to_csv(args.output_filtered, sep="\t", index=False)
```

**Line `pd.concat(rows_all, ignore_index=True)`:**
`rows_all` is a list of DataFrames — one per chunk. `pd.concat` stacks them
all into one big DataFrame. `ignore_index=True` renumbers the rows from 0
instead of keeping each chunk's internal row numbers.

**`if rows_all else pd.DataFrame(columns=...)`:**
If `rows_all` is empty (no pairs were computed at all, e.g. no known viruses),
create an empty DataFrame with the correct column names instead of crashing.

**Line `df_all.to_csv(args.output_all, sep="\t", index=False)`:**
Writes the DataFrame to a tab-separated file. `index=False` means don't write
the row numbers as a column.

### Final print (lines 376–377)

```python
print(f"Done. {len(df_all)} pairs total, {len(df_filt)} passed filters.",
      file=sys.stderr, flush=True)
```

Prints the final summary to stderr. `len(df_all)` is the total number of rows
(pairs) in the output. `flush=True` forces the message to be written immediately
(important in HPC environments where output can be buffered and appear late).

---

## The entry point (lines 380–381)

```python
if __name__ == "__main__":
    main()
```

**These two lines:**
In Python, when you run a script directly (`python script.py`), the special
variable `__name__` is set to `"__main__"`. This `if` block ensures `main()`
is only called when the script is run directly — not when it is imported as a
module by another script. This is standard Python practice for all scripts.
