import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument("--input",     required=True)
parser.add_argument("--output",    required=True)
parser.add_argument("--threshold", required=True, type=float)
args = parser.parse_args()

df  = pd.read_csv(args.input, dtype=str)
df["prob_float"] = df["prob"].astype(float)
out = df[df["prob_float"] >= args.threshold].drop(columns="prob_float")
out.to_csv(args.output, index=False)
print(f"Filtered: {len(df)} -> {len(out)} rows (prob >= {args.threshold})", file=sys.stderr)
