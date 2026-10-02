#!/bin/bash
# One-off cleanup for result/09_RDRP_structure/1_alphafold3_prediction/
# 3_data_pipeline/ batches that were left in an ambiguous state (job killed,
# node died, etc. -- no .done marker was ever written for these since the
# rule was updated to only touch it after verifying completeness).
#
# For every {source}/{batch} directory found:
#   - count input JSONs in the matching 2_batch_input/{source}/{batch}/
#   - count {name}_data.json outputs actually present in 3_data_pipeline/{source}/{batch}/
#   - if they match: touch .done (Snakemake will treat it as already built)
#   - if they don't match: remove the whole batch output directory (so
#     Snakemake reruns alphafold3_data_pipeline for it cleanly, rather than
#     resuming into a directory with some but not all data-pipeline outputs)
#
# Run from the repo root: bash scripts/09_RDRP_structure/reconcile_data_pipeline_batches.sh

set -uo pipefail

PREDDIR="result/09_RDRP_structure/1_alphafold3_prediction"
IN_ROOT="$PREDDIR/2_batch_input"
OUT_ROOT="$PREDDIR/3_data_pipeline"

if [ ! -d "$OUT_ROOT" ]; then
    echo "ERROR: $OUT_ROOT not found -- run this from the repo root" >&2
    exit 1
fi

n_complete=0
n_removed=0
n_skipped=0

for source_dir in "$OUT_ROOT"/*/; do
    source=$(basename "$source_dir")
    for batch_dir in "$source_dir"*/; do
        batch=$(basename "$batch_dir")
        [ "$batch" = "*" ] && continue

        in_dir="$IN_ROOT/$source/$batch"
        marker="$batch_dir/.done"

        if [ -e "$marker" ]; then
            n_skipped=$((n_skipped + 1))
            continue
        fi

        if [ ! -d "$in_dir" ]; then
            echo "WARNING: no matching input dir $in_dir for $source/$batch -- skipping (not touching or removing)" >&2
            n_skipped=$((n_skipped + 1))
            continue
        fi

        n_in=$(find "$in_dir" -maxdepth 1 -name '*.json' 2>/dev/null | wc -l)
        n_out=$(find "$batch_dir" -mindepth 2 -maxdepth 2 -name '*_data.json' 2>/dev/null | wc -l)

        if [ "$n_in" -gt 0 ] && [ "$n_in" -eq "$n_out" ]; then
            touch "$marker"
            echo "COMPLETE  $source/$batch  ($n_out/$n_in) -- marked done"
            n_complete=$((n_complete + 1))
        else
            echo "INCOMPLETE $source/$batch  ($n_out/$n_in) -- removing"
            rm -rf "$batch_dir"
            n_removed=$((n_removed + 1))
        fi
    done
done

echo ""
echo "Done. $n_complete marked complete, $n_removed incomplete batch dirs removed, $n_skipped already-done or unmatched skipped."
echo "Re-run snakemake to pick up the removed batches for rerun."
