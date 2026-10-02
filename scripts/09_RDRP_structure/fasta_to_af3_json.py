#!/usr/bin/env python3
"""
Convert a protein FASTA into one AlphaFold3 input JSON per sequence.

AF3's run_alphafold.py takes one job (one structure prediction) per JSON
file, so a multi-sequence FASTA is split into --outdir/{safe_id}.json, each:
    {
      "name": "{safe_id}",
      "modelSeeds": [1],
      "sequences": [{"protein": {"id": "A", "sequence": "..."}}],
      "dialect": "alphafold3",
      "version": 1
    }

Headers are sanitized to filesystem-safe job names (first whitespace-
delimited token, non-alnum/./-/_ replaced with "_").
"""

import argparse
import json
import os
import re
import sys


def read_fasta(path):
    header = None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_lines)
                header = line[1:]
                seq_lines = []
            elif header is not None:
                seq_lines.append(line)
        if header is not None:
            yield header, "".join(seq_lines)


def safe_name(header):
    seq_id = header.split()[0]
    return re.sub(r"[^A-Za-z0-9._-]", "_", seq_id)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    n = 0
    seen = set()
    for header, seq in read_fasta(args.fasta):
        seq = seq.replace("*", "").strip()
        if not seq:
            print(f"  skipping empty sequence: {header}", file=sys.stderr)
            continue
        name = safe_name(header)
        if name in seen:
            print(f"  WARNING: duplicate job name '{name}', overwriting", file=sys.stderr)
        seen.add(name)

        job = {
            "name": name,
            "modelSeeds": [args.seed],
            "sequences": [
                {"protein": {"id": "A", "sequence": seq}}
            ],
            "dialect": "alphafold3",
            "version": 1,
        }
        out_path = os.path.join(args.outdir, f"{name}.json")
        with open(out_path, "w") as out_fh:
            json.dump(job, out_fh, indent=2)
        n += 1

    print(f"Wrote {n} AF3 job JSONs to {args.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
