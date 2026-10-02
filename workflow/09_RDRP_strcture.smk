# AlphaFold3 structure prediction for full-length RdRp proteins, across
# three sequence sources:
#   new     -- your assembled-contig RdRp proteins (03's clustered set)
#              result/03_RDRP_identification/20_cluster/combined_full_c90.faa
#   esvirtu -- ESvirtu-reference RdRp proteins (03b's combined set)
#              result/03b_esvirtu_rdrp_identification/11_final/combined_full.faa
#   RVMT_ref -- RVMT top-hit reference sequences relevant to this dataset
#              (03c's Phylum-level rdrp_full.faa, concatenated -- NOT the
#              full ~77K-sequence RVMT database, just the sequences your
#              queries actually matched)
#              result/03c_ref_palm_annot/3_combined/combined_full.faa
#
# Uses the AF3 Singularity image + model weights + databases installed by
# workflow/00_alphafold3.smk (run that first). In-image script path
# confirmed via `find / -iname run_alphafold.py` inside the cford38/
# alphafold3 .sif -- it is /app/alphafold/run_alphafold.py (NOT
# /alphafold3/run_alphafold.py as the af3-on-hpc README's ${AF3_CODE_DIR}
# variable might suggest).
#
# Batching: AF3's run_alphafold.py supports --input_dir <dir> to process
# every JSON file in a directory within ONE invocation (confirmed from
# docs/input.md), rather than one job per --json_path call. Both the CPU
# data-pipeline stage and the GPU inference stage are batched at
# BATCH_SIZE sequences per SLURM job -- the CPU stage to amortize
# jackhmmer's fixed per-job startup cost, the GPU stage because DelftBlue's
# gpu-a100 partition has a wall-time limit and a single job running many
# sequences back-to-back uses that allocation more efficiently than one
# job per sequence.
#
# Job discovery: batch membership isn't known until the source FASTA has
# been split into per-sequence JSONs and grouped. This file uses the same
# plain-rule + list-file pattern as 03c_ref_palm_annot.smk /
# 101b_RDRP_phylum_tree.smk (a Snakemake checkpoint was found not to
# reliably re-trigger downstream job discovery through rule all's
# dependency chain on this Snakemake version -- see conversation notes).
# Run once to produce the per-source batch list, then again to build
# everything: see snake_run.sh for the exact two-step commands.

import os
import math

configfile: "config/config.yaml"

_AF3 = config["alphafold3"]
_OUTDIR = "result/09_RDRP_structure"
# Everything downstream of staging (JSON conversion, batching, data-pipeline,
# inference) lives under this subfolder; only 0_input stays directly under
# _OUTDIR.
_PREDDIR = _OUTDIR + "/1_alphafold3_prediction"
_BATCH_SIZE = _AF3["batch_size"]

_SOURCE_FASTA = {
    "new":      "result/03_RDRP_identification/20_cluster/combined_full_c90.faa",
    "esvirtu":  "result/03b_esvirtu_rdrp_identification/11_final/combined_full.faa",
    "RVMT_ref": "result/03c_ref_palm_annot/3_combined/combined_full.faa",
}
_SOURCES = list(_SOURCE_FASTA.keys())


def _batches_for_source(source):
    batches_file = f"{_PREDDIR}/2_batches/{source}_batches.txt"
    if not os.path.exists(batches_file):
        return []
    with open(batches_file) as fh:
        return [line.strip() for line in fh if line.strip()]


def all_structures(wildcards):
    files = []
    for source in _SOURCES:
        for batch in _batches_for_source(source):
            files.append(f"{_PREDDIR}/4_structures/{source}/{batch}/.done")
    return files


rule all:
    input:
        expand(_PREDDIR + "/2_batches/{source}_batches.txt", source=_SOURCES),
        all_structures,


# ── Step 0: stage a self-contained snapshot of each source FASTA ─────────────
# Plain copy, no conversion -- keeps this workflow's inputs independent of
# whatever else touches the upstream 03/03b/03c result directories.

rule stage_input:
    input:
        fasta = lambda wc: _SOURCE_FASTA[wc.source],
    output:
        fasta = _OUTDIR + "/0_input/{source}/source.faa",
    wildcard_constraints:
        source = "new|esvirtu|RVMT_ref",
    log:
        err = "log/09_RDRP_structure/0_input/{source}.err",
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p $(dirname {output.fasta}) $(dirname {log.err})
        cp {input.fasta} {output.fasta} 2> {log.err}
        """


# ── Step 1: convert each staged source FASTA into one AF3 input JSON per seq ─

rule fasta_to_json:
    input:
        fasta  = _OUTDIR + "/0_input/{source}/source.faa",
        script = "scripts/09_RDRP_structure/fasta_to_af3_json.py",
    output:
        outdir = directory(_PREDDIR + "/1_json/{source}"),
    wildcard_constraints:
        source = "new|esvirtu|RVMT_ref",
    log:
        out = "log/09_RDRP_structure/1_alphafold3_prediction/1_json/{source}.log",
        err = "log/09_RDRP_structure/1_alphafold3_prediction/1_json/{source}.err",
    conda:
        "../envs/python.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p {output.outdir} $(dirname {log.out})
        python {input.script} \
            --fasta  {input.fasta} \
            --outdir {output.outdir} \
            > {log.out} 2> {log.err}
        """


# ── Step 1b: group the per-sequence JSONs into batches of BATCH_SIZE ─────────
# Each batch is its own subfolder (2_batch_input/{source}/batch_NNN/) holding
# up to BATCH_SIZE JSON files, symlinked from 1_json/{source}/ -- this is
# what --input_dir points at, one call per batch instead of one call per
# sequence.

rule make_batches:
    input:
        indir = _PREDDIR + "/1_json/{source}",
    output:
        marker = touch(_PREDDIR + "/2_batch_input/{source}/.done"),
    wildcard_constraints:
        source = "new|esvirtu|RVMT_ref",
    params:
        batch_size = _BATCH_SIZE,
        outdir     = _PREDDIR + "/2_batch_input/{source}",
    log:
        err = "log/09_RDRP_structure/1_alphafold3_prediction/2_batch_input/{source}.err",
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log.err})
        : > {log.err}
        i=0
        batch=0
        for f in {input.indir}/*.json; do
            [ -e "$f" ] || continue
            if [ $((i % {params.batch_size})) -eq 0 ]; then
                batch_dir=$(printf '{params.outdir}/batch_%03d' "$batch")
                mkdir -p "$batch_dir"
                batch=$((batch + 1))
            fi
            ln -sf "$(realpath "$f")" "$batch_dir/"
            i=$((i + 1))
        done
        echo "Grouped $i JSON files into $batch batches of up to {params.batch_size}" >> {log.err}
        """


# ── Step 1c: list batch names present per source (drives job discovery) ─────

rule list_batches:
    input:
        marker = _PREDDIR + "/2_batch_input/{source}/.done",
    output:
        txt = _PREDDIR + "/2_batches/{source}_batches.txt",
    wildcard_constraints:
        source = "new|esvirtu|RVMT_ref",
    params:
        indir = _PREDDIR + "/2_batch_input/{source}",
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p $(dirname {output.txt})
        ls {params.indir} | grep '^batch_' | sort > {output.txt}
        """


# ── Step 3: CPU-only genetic/template search (--run_data_pipeline) ───────────
# jackhmmer/hmmsearch against UniRef90/UniProt/MGnify/BFD, no GPU needed --
# run on regular compute nodes. --input_dir processes every JSON in the
# batch folder in one invocation; AF3 writes one {name}_data.json per
# sequence into --output_dir.
#
# --jackhmmer_n_cpu / --nhmmer_n_cpu each independently default to
# min(cpu_count, 8) inside run_alphafold.py -- silently capped at 8
# regardless of the actual SLURM allocation, and not tied to this rule's
# `threads`/`cpus_per_task`. Pinned explicitly to data_pipeline.threads so
# raising that config value actually uses the extra CPUs instead of leaving
# them idle at the 8-CPU default. NOTE: --hmmsearch_n_cpu does NOT exist in
# this AF3 version (v3.0.0, the .sif's pinned version) -- passing it fails
# with "Unknown command line flag 'hmmsearch_n_cpu'" (confirmed via a real
# run's .err log); only jackhmmer/nhmmer expose a CPU-count flag here.

rule alphafold3_data_pipeline:
    input:
        batch_dir = _PREDDIR + "/2_batch_input/{source}/{batch}",
        sif       = _AF3["sif"],
        db_marker = _AF3["db_marker"],
    output:
        marker = _PREDDIR + "/3_data_pipeline/{source}/{batch}/.done",
    wildcard_constraints:
        source = "new|esvirtu|RVMT_ref",
    params:
        db_dir  = _AF3["db_dir"],
        outdir  = _PREDDIR + "/3_data_pipeline/{source}/{batch}",
        n_cpu   = _AF3["data_pipeline"]["threads"],
    log:
        log = "log/09_RDRP_structure/1_alphafold3_prediction/3_data_pipeline/{source}_{batch}.log",
        err = "log/09_RDRP_structure/1_alphafold3_prediction/3_data_pipeline/{source}_{batch}.err",
    threads: _AF3["data_pipeline"]["threads"]
    resources:
        mem_mb_per_cpu  = _AF3["data_pipeline"]["memory"],
        runtime         = _AF3["data_pipeline"]["runtime"],
        cpus_per_task   = _AF3["data_pipeline"]["threads"],
        slurm_partition = _AF3["data_pipeline"]["partition"],
        slurm_account   = _AF3["data_pipeline"]["account"],
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log.log})
        singularity exec {input.sif} \
            python /app/alphafold/run_alphafold.py \
                --input_dir={input.batch_dir} \
                --db_dir={params.db_dir} \
                --output_dir={params.outdir} \
                --run_data_pipeline=true \
                --run_inference=false \
                --jackhmmer_n_cpu={params.n_cpu} \
                --nhmmer_n_cpu={params.n_cpu} \
                > {log.log} 2> {log.err}

        # AF3's --input_dir call can exit 0 while silently skipping a
        # sequence it failed on (rather than aborting the whole batch) --
        # so a clean exit code alone doesn't prove every job in this batch
        # actually produced output. Only mark .done if every input JSON
        # has a matching {{name}}_data.json in the output dir.
        n_in=$(ls {input.batch_dir}/*.json 2>/dev/null | wc -l)
        n_out=$(find {params.outdir} -mindepth 2 -maxdepth 2 -name '*_data.json' 2>/dev/null | wc -l)
        if [ "$n_in" -ne "$n_out" ]; then
            echo "ERROR: expected $n_in data-pipeline outputs, found $n_out -- not marking done" >> {log.err}
            exit 1
        fi
        touch {output.marker}
        """


# ── Step 4: GPU-only structure inference (--run_inference) ───────────────────
# Consumes the data-pipeline output from Step 3 (MSAs/templates already
# computed for every sequence in this batch) -- --input_dir processes every
# {name}_data.json in that folder within one GPU allocation, batching
# BATCH_SIZE inferences per SLURM job to use gpu-a100's wall-time limit
# efficiently instead of one short job per sequence.
#
# XLA env vars per docs/performance.md's A100 80GB recommendations:
#   XLA_PYTHON_CLIENT_PREALLOCATE / XLA_CLIENT_MEM_FRACTION size the GPU
#   memory pool; XLA_FLAGS disables Triton GEMM, which the docs note causes
#   a large compile-time regression specifically on A100.
#
# No JAX compilation cache / warm-up step -- a real 10-job test run showed
# inference itself only takes ~1-2 minutes per sequence, dwarfed by the CPU
# data-pipeline stage (jackhmmer alone routinely takes 200-1400+ seconds
# per job), so caching/pre-warming compiled kernels wasn't worth the added
# complexity or a dedicated GPU allocation.

rule alphafold3_inference:
    input:
        data_marker = _PREDDIR + "/3_data_pipeline/{source}/{batch}/.done",
        sif         = _AF3["sif"],
        weights     = _AF3["model_dir"] + "/af3.bin.zst",
    output:
        marker = _PREDDIR + "/4_structures/{source}/{batch}/.done",
    wildcard_constraints:
        source = "new|esvirtu|RVMT_ref",
    params:
        model_dir = _AF3["model_dir"],
        indir     = _PREDDIR + "/3_data_pipeline/{source}/{batch}",
        outdir    = _PREDDIR + "/4_structures/{source}/{batch}",
        flat_in   = _PREDDIR + "/4_structures/{source}/{batch}/.flat_input",
    log:
        log = "log/09_RDRP_structure/1_alphafold3_prediction/4_structures/{source}_{batch}.log",
        err = "log/09_RDRP_structure/1_alphafold3_prediction/4_structures/{source}_{batch}.err",
    threads: _AF3["threads"]
    resources:
        mem_mb_per_cpu  = _AF3["memory"],
        runtime         = _AF3["runtime"],
        cpus_per_task   = _AF3["threads"],
        slurm_partition = _AF3["partition"],
        slurm_account   = _AF3["account"],
        # `gpus_per_task` is NOT a resource key the snakemake-executor-
        # plugin-slurm recognizes -- it was silently dropped, so past runs
        # submitted with no GPU request at all (confirmed by a real job's
        # .err: CUDA_ERROR_NO_DEVICE / "No visible GPU devices" even though
        # it landed on a gpu-a100 node). DelftBlue's sbatch separately
        # rejects a bare `--gpus` flag without `--ntasks` ("Use '--gpus-
        # per-task' instead" -- confirmed via test/test_af3_10jobs.sh,
        # which works with --ntasks=1 --gpus-per-task=1). Route both
        # through slurm_extra, the plugin's supported passthrough for
        # flags it doesn't natively model.
        slurm_extra     = f"--ntasks=1 --gpus-per-task={_AF3['gpus']}",
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log.log})

        # AF3's --input_dir only looks at the immediate top level of the
        # given directory, NOT recursively into subfolders (confirmed
        # empirically: pointing it at 3_data_pipeline/{{source}}/{{batch}}
        # directly -- which holds {{name}}/{{name}}_data.json one level
        # deep -- reported "Processing 0 fold inputs" despite 25 real
        # *_data.json files being present). Flatten by symlinking each
        # {{name}}/{{name}}_data.json into a scratch dir at the top level
        # before pointing --input_dir at that instead.
        #
        # Resumable, same pattern as alphafold3_data_pipeline: skip
        # sequences that already have a {{name}}/{{name}}_model.cif from a
        # prior (killed) attempt.
        rm -rf {params.flat_in}
        mkdir -p {params.flat_in}
        n_skipped=0
        for f in {params.indir}/*/*_data.json; do
            [ -e "$f" ] || continue
            name=$(basename "$f" _data.json)
            if [ -e "{params.outdir}/$name/${{name}}_model.cif" ]; then
                n_skipped=$((n_skipped + 1))
            else
                ln -sf "$(realpath "$f")" "{params.flat_in}/"
            fi
        done
        n_remaining=$(ls {params.flat_in}/*.json 2>/dev/null | wc -l)
        echo "$n_skipped already complete from a prior attempt, $n_remaining remaining to run" > {log.log}

        export CUDA_VISIBLE_DEVICES={_AF3[gpu_devices]}
        export XLA_PYTHON_CLIENT_PREALLOCATE=true
        export XLA_CLIENT_MEM_FRACTION=0.95
        export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"

        if [ "$n_remaining" -gt 0 ]; then
            singularity exec --nv {input.sif} \
                python /app/alphafold/run_alphafold.py \
                    --input_dir={params.flat_in} \
                    --model_dir={params.model_dir} \
                    --output_dir={params.outdir} \
                    --run_data_pipeline=false \
                    --run_inference=true \
                    >> {log.log} 2> {log.err}
        fi
        rm -rf {params.flat_in}

        # Same guard as alphafold3_data_pipeline: a clean exit code doesn't
        # prove every sequence in this batch actually got a structure --
        # only mark .done if every input {{name}}_data.json (across all
        # attempts) has a matching {{name}}_model.cif in the output dir.
        n_in=$(find {params.indir} -mindepth 2 -maxdepth 2 -name '*_data.json' 2>/dev/null | wc -l)
        n_out=$(find {params.outdir} -mindepth 2 -maxdepth 2 -name '*_model.cif' 2>/dev/null | wc -l)
        if [ "$n_in" -ne "$n_out" ]; then
            echo "ERROR: expected $n_in structures, found $n_out -- not marking done" >> {log.err}
            exit 1
        fi
        touch {output.marker}
        """
