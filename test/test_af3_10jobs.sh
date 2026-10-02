#!/bin/bash
#SBATCH --job-name=af3_test10
# gpu-v100 (Tesla V100S, CC 7.0) is too old for AF3 (needs CC >=8.0) --
# gpu-a100 (Tesla A100 80GB, CC 8.0) is required.
#SBATCH --partition=gpu-a100
#SBATCH --account=research-ceg-wm
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8000
#SBATCH --time=06:00:00
#SBATCH --output=log/09_RDRP_structure/test_af3_10jobs_%j.log
#SBATCH --error=log/09_RDRP_structure/test_af3_10jobs_%j.err

# Quick smoke test: run AF3 on the first 10 JSON jobs found (across whichever
# sources under result/09_RDRP_structure/1_json/ already exist), skipping
# Snakemake entirely. Run 09_RDRP_strcture.smk's stage_input + fasta_to_json
# rules (or the whole pass-1 snakemake command in snake_run.sh) first so
# 1_json/{source}/*.json actually exist before submitting this.
#
# Submit with: sbatch test/test_af3_10jobs.sh

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

SIF="tools/alphafold3/alphafold3.sif"
MODEL_DIR="database/alphafold3/models"
DB_DIR="database/alphafold3/public_databases"
JSON_ROOT="result/09_RDRP_structure/1_json"
OUT_ROOT="result/09_RDRP_structure/test_10job_run"

mkdir -p "$OUT_ROOT" log/09_RDRP_structure

# Collect up to 10 JSON files across all sources (new/esvirtu/RVMT_ref),
# deterministic order (sorted by path).
mapfile -t JOBS < <(find "$JSON_ROOT" -name '*.json' | sort | head -n 10)

if [ "${#JOBS[@]}" -eq 0 ]; then
    echo "No JSON files found under $JSON_ROOT -- run the fasta_to_json step first." >&2
    exit 1
fi

echo "Testing ${#JOBS[@]} job(s):"
printf '  %s\n' "${JOBS[@]}"

export CUDA_VISIBLE_DEVICES=0

for json_file in "${JOBS[@]}"; do
    source_name=$(basename "$(dirname "$json_file")")
    job_name=$(basename "$json_file" .json)
    outdir="$OUT_ROOT/$source_name"
    mkdir -p "$outdir"

    echo "=== $source_name / $job_name ==="
    singularity exec --nv "$SIF" \
        python /app/alphafold/run_alphafold.py \
            --json_path="$json_file" \
            --model_dir="$MODEL_DIR" \
            --db_dir="$DB_DIR" \
            --output_dir="$outdir"
done

echo "Done. Structures (if successful) under $OUT_ROOT/{source}/{job}/{job}_model.cif"
