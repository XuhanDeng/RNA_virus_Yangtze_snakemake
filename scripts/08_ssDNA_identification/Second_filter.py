import argparse
import os
from typing import Tuple

import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Second filter: VS2 hallmark/full + GeNomad non-provirus + CheckV viral_genes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i1", "--checkv", help="quality_summary.tsv")
    parser.add_argument("-i2", "--vs2", help="virsorter_result;final-viral-score.tsv")
    parser.add_argument("-i3", "--genomad", help="virsorter_result;_virus_summary.tsv")
    parser.add_argument("-o", "--output", help="the output dir to store files")
    return parser.parse_args()


def load_inputs(checkv_path: str, vs2_path: str, genomad_path: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load CheckV, VirSorter2, and GeNomad input tables."""
    checkv_df = pd.read_csv(checkv_path, sep="\t")
    vs2_df = pd.read_csv(vs2_path, sep="\t")
    genomad_df = pd.read_csv(genomad_path, sep="\t")
    return checkv_df, vs2_df, genomad_df


def filter_vs2(vs2_df: pd.DataFrame) -> pd.DataFrame:
    """Filter VirSorter2 results: hallmark > 0 and full sequences."""
    vs2_df = vs2_df.copy()
    vs2_df["name"] = vs2_df["seqname"].str.split("|").str[0]
    vs2_df["full/partial"] = vs2_df["seqname"].str.split("|").str[2]
    vs2_filter = vs2_df[(vs2_df["hallmark"] > 0) & (vs2_df["full/partial"] == "full")]
    return vs2_filter.loc[:, ["name"]]


def filter_genomad(genomad_df: pd.DataFrame) -> pd.DataFrame:
    """Filter GeNomad results: exclude provirus and require hallmark genes."""
    genomad_df = genomad_df.copy()
    genomad_df["name"] = genomad_df["seq_name"]
    genomad_filter = genomad_df[(genomad_df["topology"] != "Provirus") & (genomad_df["n_hallmarks"] > 0)]
    return genomad_filter.loc[:, ["name"]]


def filter_checkv(checkv_df: pd.DataFrame) -> pd.DataFrame:
    """Filter CheckV results: viral_genes > 0."""
    checkv_df = checkv_df.copy()
    checkv_df["name"] = checkv_df["contig_id"]
    checkv_filter = checkv_df[checkv_df["viral_genes"] > 0]
    return checkv_filter.loc[:, ["name"]]


def merge_lists(vs2_list: pd.DataFrame, genomad_list: pd.DataFrame, checkv_list: pd.DataFrame) -> pd.DataFrame:
    """Merge name lists and remove duplicates."""
    merged = pd.DataFrame()
    merged["name"] = pd.concat(
        [vs2_list.iloc[:, 0], genomad_list.iloc[:, 0], checkv_list.iloc[:, 0]],
        ignore_index=True,
    )
    return merged.drop_duplicates("name")


def build_source_table(
    vs2_list: pd.DataFrame, genomad_list: pd.DataFrame, checkv_list: pd.DataFrame
) -> pd.DataFrame:
    """Build a table indicating which tool(s) identified each sequence."""
    names = pd.DataFrame(
        {
            "name": pd.concat(
                [vs2_list.iloc[:, 0], genomad_list.iloc[:, 0], checkv_list.iloc[:, 0]],
                ignore_index=True,
            ).drop_duplicates()
        }
    )
    vs2_flag = vs2_list.drop_duplicates("name").assign(vs2=True)
    genomad_flag = genomad_list.drop_duplicates("name").assign(genomad=True)
    checkv_flag = checkv_list.drop_duplicates("name").assign(checkv=True)

    merged = names.merge(vs2_flag, on="name", how="left")
    merged = merged.merge(genomad_flag, on="name", how="left")
    merged = merged.merge(checkv_flag, on="name", how="left")
    merged[["vs2", "genomad", "checkv"]] = merged[["vs2", "genomad", "checkv"]].fillna(False)
    merged["identified_by"] = (
        merged[["vs2", "genomad", "checkv"]]
        .apply(lambda row: "|".join([k for k, v in row.items() if v]), axis=1)
        .replace("", "none")
    )
    return merged


def write_output(
    names: pd.DataFrame,
    base_name: str,
    output_dir: str,
    vs2_filtered: pd.DataFrame,
    genomad_filtered: pd.DataFrame,
    checkv_filtered: pd.DataFrame,
    source_table: pd.DataFrame,
) -> None:
    """Write filtered list and summary tables to output directory."""
    os.makedirs(output_dir, exist_ok=True)
    names.to_csv(os.path.join(output_dir, base_name + "_checkv_extract.txt"), index=False, header=None)
    vs2_filtered.to_csv(os.path.join(output_dir, base_name + "_vs2_filtered.tsv"), sep="\t", index=False)
    genomad_filtered.to_csv(os.path.join(output_dir, base_name + "_genomad_filtered.tsv"), sep="\t", index=False)
    checkv_filtered.to_csv(os.path.join(output_dir, base_name + "_checkv_filtered.tsv"), sep="\t", index=False)
    source_table.to_csv(os.path.join(output_dir, base_name + "_source_table.tsv"), sep="\t", index=False)


def main() -> None:
    """Run second-pass filtering of viral contigs."""
    args = parse_args()
    checkv_df, vs2_df, genomad_df = load_inputs(args.checkv, args.vs2, args.genomad)
    vs2_list = filter_vs2(vs2_df)
    genomad_list = filter_genomad(genomad_df)
    checkv_list = filter_checkv(checkv_df)
    merged = merge_lists(vs2_list, genomad_list, checkv_list)
    source_table = build_source_table(vs2_list, genomad_list, checkv_list)
    base_name = os.path.basename(args.genomad).split("_scaffolds_")[0]
    write_output(
        merged,
        base_name,
        args.output,
        vs2_list,
        genomad_list,
        checkv_list,
        source_table,
    )


if __name__ == "__main__":
    main()
