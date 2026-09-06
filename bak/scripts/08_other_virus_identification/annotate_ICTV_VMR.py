#!/usr/bin/env python3
"""
Annotate GeNomad virus_filter_summary with ICTV VMR taxonomy and genome composition.

GeNomad 'taxonomy' column is a semicolon-separated lineage string:
  "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;Ackermannviridae"
Empty fields (;;) are skipped. Taxa are matched right-to-left (most specific first)
against VMR Family/Order/Class/Phylum/Kingdom/Realm columns.

When a taxon maps to multiple VMR rows with different Genome composition or Host Source
values (e.g. Circoviridae has ssDNA(+) and ssDNA(-) members), all unique values are
joined with "; " so no information is lost.

When matched at a higher rank (e.g. Class), more specific columns (Order, Family) are
left empty to avoid reporting a spurious lower-rank value.
"""

import argparse
import re
import sys

import pandas as pd


# Rank columns from most specific to most general
RANK_COLS = ["Family", "Order", "Class", "Phylum", "Kingdom", "Realm"]

VMR_OUT_COLS = [
    "Realm", "Kingdom", "Phylum", "Class", "Order", "Family",
    "Genome composition", "Host Source"
]

def unify_genome_composition(value: str) -> str:
    """Strip strand polarity suffixes and deduplicate base types.

    'ssDNA; ssDNA(+); ssDNA(+/-)' -> 'ssDNA'
    'ssDNA; ssRNA'                 -> 'ssDNA; ssRNA'
    ''                             -> ''
    """
    if not value or not value.strip():
        return ""
    parts = [v.strip() for v in value.split(";") if v.strip()]
    # Remove trailing (...) polarity annotation, case-insensitive
    bases = [re.sub(r"\s*\([^)]*\)\s*$", "", p).strip() for p in parts]
    seen = []
    for b in bases:
        if b and b not in seen:
            seen.append(b)
    return "; ".join(seen)


def load_vmr(xlsx: str, sheet: str) -> pd.DataFrame:
    df = pd.read_excel(xlsx, sheet_name=sheet, dtype=str)
    df.columns = df.columns.str.strip()
    keep = [c for c in RANK_COLS + ["Genome composition", "Host Source"] if c in df.columns]
    df = df[keep].fillna("")
    print(f"VMR loaded: {len(df)} rows", file=sys.stderr)
    return df


def build_lookup(vmr: pd.DataFrame) -> dict:
    """
    Map taxon_name -> (matched_rank_col, aggregated_dict).

    Taxonomy columns (Realm..Family) take the first non-empty value since they
    are consistent within a taxon. Genome composition and Host Source are
    aggregated as sorted unique values joined with '; ' to handle mixed families.
    """
    lookup = {}
    for col in RANK_COLS:
        if col not in vmr.columns:
            continue
        groups = vmr[vmr[col] != ""].groupby(col)
        for taxon, grp in groups:
            taxon = taxon.strip()
            if not taxon or taxon in lookup:
                continue
            matched_idx = RANK_COLS.index(col)
            tax_dict = {}
            # taxonomy columns: consistent within group, take first non-empty
            for i, rc in enumerate(RANK_COLS):
                if i >= matched_idx and rc in grp.columns:
                    vals = grp[rc].replace("", pd.NA).dropna().unique()
                    tax_dict[rc] = str(vals[0]).strip() if len(vals) > 0 else ""
                else:
                    tax_dict[rc] = ""  # below match rank — leave empty
            # aggregate columns: join all unique non-empty values
            for agg_col in ("Genome composition", "Host Source"):
                if agg_col in grp.columns:
                    unique_vals = sorted(
                        v.strip() for v in grp[agg_col].unique()
                        if str(v).strip() and str(v).strip() != "nan"
                    )
                    tax_dict[agg_col] = "; ".join(unique_vals)
                else:
                    tax_dict[agg_col] = ""
            lookup[taxon] = (col, tax_dict)
    return lookup


def annotate_row(taxonomy: str, lookup: dict) -> dict:
    empty = {c: "" for c in VMR_OUT_COLS}
    if not taxonomy or str(taxonomy).strip() in ("nan", "Unclassified"):
        return empty
    taxa = [t.strip() for t in str(taxonomy).split(";") if t.strip()]
    for taxon in reversed(taxa):
        if taxon in lookup:
            _, tax_dict = lookup[taxon]
            return {c: tax_dict.get(c, "") for c in VMR_OUT_COLS}
    return empty


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--vmr",     required=True)
    parser.add_argument("--sheet",   required=True)
    parser.add_argument("--output",  required=True)
    args = parser.parse_args()

    vmr    = load_vmr(args.vmr, args.sheet)
    lookup = build_lookup(vmr)
    print(f"Lookup: {len(lookup)} unique taxon names", file=sys.stderr)

    df = pd.read_csv(args.summary, sep="\t", dtype=str)

    if df.empty or "taxonomy" not in df.columns:
        for c in VMR_OUT_COLS + ["Unified Genome Composition"]:
            df[c] = ""
        df.to_csv(args.output, sep="\t", index=False)
        print("Empty input or missing taxonomy column", file=sys.stderr)
        return

    annot_df = pd.DataFrame(list(df["taxonomy"].apply(lambda t: annotate_row(t, lookup))))
    annot_df["Unified Genome Composition"] = annot_df["Genome composition"].apply(unify_genome_composition)
    df = pd.concat([df, annot_df], axis=1)
    df.to_csv(args.output, sep="\t", index=False)

    matched = (annot_df["Realm"] != "").sum()
    print(f"Annotated {matched}/{len(df)} rows -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
