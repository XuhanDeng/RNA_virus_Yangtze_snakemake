#!/usr/bin/env python3
"""Filter viral contigs using GeNomad + CheckV (union logic, no VirSorter2)."""

import argparse
import os

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter: GeNomad non-provirus+hallmark OR CheckV viral_genes>0"
    )
    parser.add_argument("--checkv",  required=True, help="CheckV quality_summary.tsv")
    parser.add_argument("--genomad", required=True, help="GeNomad *_virus_summary.tsv")
    parser.add_argument("--output",  required=True, help="Output directory")
    return parser.parse_args()


def filter_genomad(genomad_df: pd.DataFrame) -> pd.DataFrame:
    """Keep non-provirus sequences with at least one hallmark gene."""
    df = genomad_df.copy()
    df["name"] = df["seq_name"]
    filtered = df[(df["topology"] != "Provirus") & (df["n_hallmarks"] > 0)]
    return filtered[["name"]]


def filter_checkv(checkv_df: pd.DataFrame) -> pd.DataFrame:
    """Keep sequences with at least one viral gene detected."""
    df = checkv_df.copy()
    df["name"] = df["contig_id"]
    filtered = df[df["viral_genes"] > 0]
    return filtered[["name"]]


def build_source_table(genomad_list: pd.DataFrame, checkv_list: pd.DataFrame) -> pd.DataFrame:
    all_names = pd.concat(
        [genomad_list["name"], checkv_list["name"]], ignore_index=True
    ).drop_duplicates()
    names = pd.DataFrame({"name": all_names})
    genomad_flag = genomad_list.drop_duplicates("name").assign(genomad=True)
    checkv_flag  = checkv_list.drop_duplicates("name").assign(checkv=True)
    merged = names.merge(genomad_flag, on="name", how="left")
    merged = merged.merge(checkv_flag,  on="name", how="left")
    merged[["genomad", "checkv"]] = merged[["genomad", "checkv"]].fillna(False)
    merged["identified_by"] = (
        merged[["genomad", "checkv"]]
        .apply(lambda row: "|".join(k for k, v in row.items() if v), axis=1)
        .replace("", "none")
    )
    return merged


def main() -> None:
    args = parse_args()
    genomad_df = pd.read_csv(args.genomad, sep="\t")
    checkv_df  = pd.read_csv(args.checkv,  sep="\t")

    genomad_list = filter_genomad(genomad_df)
    checkv_list  = filter_checkv(checkv_df)

    merged = pd.concat(
        [genomad_list, checkv_list], ignore_index=True
    ).drop_duplicates("name")

    source_table = build_source_table(genomad_list, checkv_list)

    base_name = os.path.basename(args.genomad).split("_scaffolds_")[0]
    os.makedirs(args.output, exist_ok=True)

    merged.to_csv(
        os.path.join(args.output, base_name + "_checkv_extract.txt"),
        index=False, header=None
    )
    genomad_list.to_csv(
        os.path.join(args.output, base_name + "_genomad_filtered.tsv"),
        sep="\t", index=False
    )
    checkv_list.to_csv(
        os.path.join(args.output, base_name + "_checkv_filtered.tsv"),
        sep="\t", index=False
    )
    source_table.to_csv(
        os.path.join(args.output, base_name + "_source_table.tsv"),
        sep="\t", index=False
    )

    print(f"{base_name}: genomad={len(genomad_list)}, checkv={len(checkv_list)}, union={len(merged)}")


if __name__ == "__main__":
    main()
