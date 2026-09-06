#!/usr/bin/env python3
"""
Rename tip labels in a Newick tree to shorter, readable IDs.

Rules (applied in order):
  1. LucaProt:   ORF{N}_{sample}_{contig}:{start}:{end}  -> {sample}_{contig}
  2. RdRpCATCH / ESvirtu: {prefix}_frame=N_RdRp_X-Y      -> {prefix}
  3. Everything else: left unchanged (RT outgroup, ICTV, etc.)

A lookup TSV is written alongside the renamed tree:
  original_id <TAB> renamed_id
"""

import re
import sys


def rename_tip(tip: str) -> str:
    # LucaProt: ORF{N}_{rest}:{start}:{end}
    m = re.match(r"^ORF\d+_(.+?)(?::\d+:\d+)?$", tip)
    if m:
        return m.group(1)

    # RdRpCATCH / ESvirtu: strip _frame=... suffix
    m = re.match(r"^(.+?)_frame=.*$", tip)
    if m:
        return m.group(1)

    return tip


def process_newick(nwk_text: str) -> tuple[str, dict]:
    lookup = {}

    def replacer(match):
        original = match.group(0)
        renamed  = rename_tip(original)
        lookup[original] = renamed
        return renamed

    # Match tip labels: must start with a non-digit, non-dot character
    # This excludes branch lengths (e.g. 0.000123) from being matched
    renamed = re.sub(r"[A-Za-z_][^(),;:]*(?=[:(),;])", replacer, nwk_text)
    return renamed, lookup


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",   required=True, help="input .nwk file")
    parser.add_argument("--output",  required=True, help="output renamed .nwk file")
    parser.add_argument("--lookup",  required=True, help="output lookup TSV (original -> renamed)")
    args = parser.parse_args()

    with open(args.input) as fh:
        nwk_text = fh.read()

    renamed_nwk, lookup = process_newick(nwk_text)

    with open(args.output, "w") as fh:
        fh.write(renamed_nwk)

    with open(args.lookup, "w") as fh:
        fh.write("original_id\trenamed_id\n")
        for orig, new in sorted(lookup.items()):
            fh.write(f"{orig}\t{new}\n")

    renamed_count = sum(1 for o, r in lookup.items() if o != r)
    print(f"Tips total: {len(lookup)}, renamed: {renamed_count}", file=sys.stderr)
    print(f"Written: {args.output}", file=sys.stderr)
    print(f"Written: {args.lookup}", file=sys.stderr)


if __name__ == "__main__":
    main()
