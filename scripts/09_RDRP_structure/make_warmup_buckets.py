#!/usr/bin/env python3
"""
Generate one synthetic warm-up AF3 job per distinct token-length bucket
actually present across a set of AF3 input JSONs.

AF3 pads each input to the smallest --buckets boundary >= its token count
and JAX compiles a separate kernel per bucket (run_alphafold.py's _BUCKETS
flag). Warming up only one sequence length primes one bucket; real batches
with different lengths still pay full compile cost on first use. This
script inspects every job JSON under --json_dirs, buckets its sequence
length using the same --buckets boundaries the real runs will use, and
writes one placeholder job per bucket actually needed -- so warm-up
compiles exactly the buckets the real data will hit, no more, no less.

Note: bucketing here is approximated from raw protein sequence length,
while AF3's actual token count includes MSA/template features added
during the data pipeline -- so this is a close but not exact proxy, and a
placeholder padded to a given bucket may compile a slightly different
final token count than a real sequence assigned to the same nominal
bucket. For sequences longer than the largest --buckets boundary, AF3
creates an exact-size bucket per input rather than reusing the largest one
(per its own flag docs) -- since that exact size can't be predicted ahead
of the data pipeline, such sequences are reported but NOT pre-warmed; they
will compile fresh on first real use regardless.

Placeholder sequences are a repeating motif long enough to fill each
bucket (capped at the bucket size); AF3 only cares about token count for
compilation, not biological content.
"""

import argparse
import glob
import json
import os
import sys

DEFAULT_BUCKETS = [128, 256, 384, 512, 768, 1024, 1280, 1536, 2048, 2560,
                    3072, 3584, 4096, 4608, 5120]

_MOTIF = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWELVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQAGVWPAAVRESVPSLL"


def bucket_for_length(n, buckets):
    """Return the matching --buckets boundary, or None if n exceeds all of
    them (AF3 would create an exact-size bucket for n itself in that case --
    not predictable/pre-warmable here, see module docstring)."""
    for b in buckets:
        if n <= b:
            return b
    return None


def seq_length(job):
    total = 0
    for entry in job.get("sequences", []):
        prot = entry.get("protein")
        if prot:
            total += len(prot.get("sequence", ""))
    return total


def placeholder_sequence(bucket):
    reps = bucket // len(_MOTIF) + 1
    return (_MOTIF * reps)[:bucket]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json-dirs", nargs="+", required=True,
                         help="Directories of *.json AF3 input jobs to scan for lengths")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--buckets", default=",".join(str(b) for b in DEFAULT_BUCKETS),
                         help="Comma-separated AF3 bucket boundaries (must match --buckets passed to run_alphafold.py)")
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    buckets = sorted(int(b) for b in args.buckets.split(","))
    os.makedirs(args.outdir, exist_ok=True)

    needed_buckets = set()
    n_scanned = 0
    n_oversized = 0
    for d in args.json_dirs:
        for path in glob.glob(os.path.join(d, "*.json")):
            with open(path) as fh:
                job = json.load(fh)
            n = seq_length(job)
            if n <= 0:
                continue
            n_scanned += 1
            bucket = bucket_for_length(n, buckets)
            if bucket is None:
                n_oversized += 1
                continue
            needed_buckets.add(bucket)

    if n_oversized:
        print(f"WARNING: {n_oversized} sequence(s) exceed the largest bucket "
              f"({buckets[-1]}) -- AF3 will compile an exact-size bucket for "
              f"each on first real use; not pre-warmed here.", file=sys.stderr)

    if not needed_buckets:
        print("No sequences found to bucket -- defaulting to the smallest bucket only", file=sys.stderr)
        needed_buckets = {buckets[0]}

    for bucket in sorted(needed_buckets):
        job = {
            "name": f"warmup_{bucket}",
            "modelSeeds": [args.seed],
            "sequences": [
                {"protein": {"id": "A", "sequence": placeholder_sequence(bucket)}}
            ],
            "dialect": "alphafold3",
            "version": 1,
        }
        out_path = os.path.join(args.outdir, f"warmup_{bucket}.json")
        with open(out_path, "w") as out_fh:
            json.dump(job, out_fh, indent=2)

    print(f"Scanned {n_scanned} sequences across {len(args.json_dirs)} dir(s); "
          f"wrote {len(needed_buckets)} warm-up job(s) for buckets: {sorted(needed_buckets)}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
