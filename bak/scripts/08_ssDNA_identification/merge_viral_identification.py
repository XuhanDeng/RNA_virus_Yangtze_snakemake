import argparse
import os
from typing import Tuple
import pandas as pd


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="virsoter2,max_score>=0.9;genomad != provirus",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i2", "--vs2", help="virsorter_result;final-viral-score.tsv")
    parser.add_argument("-i3", "--genomad", help="virsorter_result;_virus_summary.tsv")
    parser.add_argument("-o", "--output", help="the out put dir you want store file")
    parser.add_argument(
        "--vs2-min-score",
        type=float,
        default=0.9,
        help="Minimum VirSorter2 max_score threshold",
    )
    parser.add_argument(
        "--vs2-allow-partial",
        action="store_true",
        default=False,
        help="Allow VirSorter2 entries with full/partial != full",
    )
    parser.add_argument(
        "--genomad-exclude-topology",
        default="Provirus",
        help="Exclude GeNomad entries with this topology value",
    )
    return parser.parse_args()


def load_inputs(vs2_path: str, genomad_path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load VS2 and GeNomad input tables."""
    vs2_df = pd.read_csv(vs2_path, sep="\t")
    genomad_df = pd.read_csv(genomad_path, sep="\t")
    return vs2_df, genomad_df


def filter_vs2(vs2_df: pd.DataFrame, min_score: float, allow_partial: bool) -> pd.DataFrame:
    """Filter VirSorter2 results and return name/tag table."""
    vs2_df = vs2_df.copy()
    vs2_df["tag2"] = "vs2"
    vs2_df["name"] = vs2_df["seqname"].str.split("|").str[0]
    vs2_df["full/partial"] = vs2_df["seqname"].str.split("|").str[2]
    vs2_filter = vs2_df[vs2_df["max_score"] >= min_score]
    if not allow_partial:
        vs2_filter = vs2_filter[vs2_filter["full/partial"] == "full"]
    return vs2_filter.loc[:, ["name", "tag2"]]


def filter_genomad(genomad_df: pd.DataFrame, exclude_topology: str) -> pd.DataFrame:
    """Filter GeNomad results and return name/tag table."""
    genomad_df = genomad_df.copy()
    genomad_df["tag3"] = "genomad"
    genomad_df["name"] = genomad_df["seq_name"]
    genomad_filter = genomad_df[genomad_df["topology"] != exclude_topology]
    return genomad_filter.loc[:, ["name", "tag3"]]


def merge_labels(genomad_list: pd.DataFrame, vs2_list: pd.DataFrame) -> pd.DataFrame:
    """Merge name lists and attach source labels."""
    merge_index = pd.DataFrame()
    merge_index["name"] = pd.concat(
        [genomad_list.loc[:, "name"], vs2_list.loc[:, "name"]], ignore_index=True
    )
    merge_unique = merge_index.drop_duplicates("name")
    merge_temp = pd.merge(merge_unique, genomad_list, how="left", on="name")
    return pd.merge(merge_temp, vs2_list, how="left", on="name")


def write_outputs(merged: pd.DataFrame, base_name: str, output_dir: str) -> None:
    """Write merged CSV and sequence list outputs."""
    csv_dir = os.path.join(output_dir, "csv")
    list_dir = os.path.join(output_dir, "list")
    os.makedirs(csv_dir, exist_ok=True)
    os.makedirs(list_dir, exist_ok=True)

    merged.to_csv(os.path.join(csv_dir, base_name + "_merged_results.csv"), index=False)
    merged.loc[:, "name"].to_csv(
        os.path.join(list_dir, base_name + "_merge3_list.txt"),
        index=False,
        header=None,
    )


def main() -> None:
    """Merge GeNomad and VirSorter2 results into a combined label table."""
    args = parse_args()
    vs2_df, genomad_df = load_inputs(args.vs2, args.genomad)
    vs2_list = filter_vs2(vs2_df, args.vs2_min_score, args.vs2_allow_partial)
    genomad_list = filter_genomad(genomad_df, args.genomad_exclude_topology)
    merged = merge_labels(genomad_list, vs2_list)
    base_name = os.path.basename(args.genomad).split("_scaffolds_")[0]
    write_outputs(merged, base_name, args.output)


if __name__ == "__main__":
    main()
