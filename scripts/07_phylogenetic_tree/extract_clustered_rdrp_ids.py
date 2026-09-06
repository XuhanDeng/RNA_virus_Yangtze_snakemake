#!/usr/bin/env python3
"""
Extract protein IDs from RdRpCATCH/LucaProt FAA files whose parent contig
appears in rna_virus_contig_table.tsv (i.e. passed the 99_tables filters).
"""

import re
import sys
import pandas as pd

table_path  = snakemake.input.table
merged_path = snakemake.input.merged
faa_paths   = [snakemake.input.rc_faa, snakemake.input.lp_faa]
out_path    = snakemake.output.ids


def to_contig(pid: str) -> str:
    pid = re.sub(r"^ORF\d+_", "", pid)
    pid = pid.split(":")[0]
    pid = re.sub(r"_frame=.*", "", pid)
    return pid


table  = pd.read_csv(table_path,  sep="\t", usecols=["contig_id"])
merged = pd.read_csv(merged_path, sep="\t", usecols=["contig", "source"])

clustered = set(table["contig_id"])
matched   = merged[merged["contig"].isin(clustered)]
print(f"Clustered contigs: {len(clustered)}, matched in merged: {len(matched)}", file=sys.stderr)

ids = []
for faa_path in faa_paths:
    with open(faa_path) as fh:
        for line in fh:
            if line.startswith(">"):
                pid = line[1:].split()[0]
                if to_contig(pid) in clustered:
                    ids.append(pid)

with open(out_path, "w") as out:
    for pid in ids:
        out.write(pid + "\n")

print(f"Protein IDs written: {len(ids)}", file=sys.stderr)
