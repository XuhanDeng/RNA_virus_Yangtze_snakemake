#!/bin/bash
#SBATCH --job-name=af3_infer_test
# gpu-v100 (Tesla V100S, CC 7.0) is too old for AF3 (needs CC >=8.0) --
# gpu-a100 (Tesla A100 80GB, CC 8.0) is required.
#SBATCH --partition=gpu-a100
#SBATCH --account=research-ceg-wm
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8000
#SBATCH --time=06:00:00
#SBATCH --output=log/09_RDRP_structure/test_af3_inference_batch001_%j.log
#SBATCH --error=log/09_RDRP_structure/test_af3_inference_batch001_%j.err

# GPU-inference-only smoke test: skips the data-pipeline stage entirely,
# reusing the already-completed data-pipeline outputs from
# result/09_RDRP_structure/1_alphafold3_prediction/3_data_pipeline/new/batch_001
# (confirmed finished -- every sequence in that batch has its
# {name}/{name}_data.json already on disk). Runs run_alphafold.py with
# --run_data_pipeline=false --run_inference=true.
#
# IMPORTANT: --input_dir only looks at the immediate top level of the given
# directory, NOT recursively into subfolders -- confirmed empirically: an
# earlier version of this script pointed --input_dir directly at
# batch_001/ (which holds {name}/{name}_data.json one level deep) and AF3
# reported "Processing 0 fold inputs" despite 25 real *_data.json files
# being present. This version flattens by symlinking each
# {name}/{name}_data.json into a scratch dir at the top level first, then
# points --input_dir at that flat scratch dir instead. This is the first
# real test of whether a bare symlinked *_data.json (with no sibling
# files/subdirectory alongside it) is sufficient for inference, or whether
# AF3 needs other data-pipeline output files to travel with it -- if this
# still produces 0 structures, that's the next thing to investigate.
#
# Submit with: sbatch test/test_af3_inference_batch001.sh

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

SIF="tools/alphafold3/alphafold3.sif"
MODEL_DIR="database/alphafold3/models"
IN_DIR="result/09_RDRP_structure/1_alphafold3_prediction/3_data_pipeline/new/batch_001"
OUT_DIR="result/09_RDRP_structure/test_inference_batch001"
FLAT_IN="$OUT_DIR/.flat_input"

mkdir -p "$OUT_DIR" log/09_RDRP_structure

if [ ! -d "$IN_DIR" ]; then
    echo "ERROR: $IN_DIR not found -- run the data-pipeline stage for this batch first." >&2
    exit 1
fi

n_in=$(find "$IN_DIR" -mindepth 2 -maxdepth 2 -name '*_data.json' | wc -l)
echo "Found $n_in completed data-pipeline output(s) in $IN_DIR"
if [ "$n_in" -eq 0 ]; then
    echo "ERROR: no *_data.json files found under $IN_DIR -- nothing to run inference on." >&2
    exit 1
fi

rm -rf "$FLAT_IN"
mkdir -p "$FLAT_IN"
for f in "$IN_DIR"/*/*_data.json; do
    [ -e "$f" ] || continue
    ln -sf "$(realpath "$f")" "$FLAT_IN/"
done
n_flat=$(find "$FLAT_IN" -maxdepth 1 -name '*_data.json' | wc -l)
echo "Flattened $n_flat *_data.json symlink(s) into $FLAT_IN"

export CUDA_VISIBLE_DEVICES=0
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"

singularity exec --nv "$SIF" \
    python /app/alphafold/run_alphafold.py \
        --input_dir="$FLAT_IN" \
        --model_dir="$MODEL_DIR" \
        --output_dir="$OUT_DIR" \
        --run_data_pipeline=false \
        --run_inference=true

n_out=$(find "$OUT_DIR" -mindepth 2 -maxdepth 2 -name '*_model.cif' | wc -l)
echo "Done. $n_out/$n_in structure(s) written under $OUT_DIR/{name}/{name}_model.cif"
