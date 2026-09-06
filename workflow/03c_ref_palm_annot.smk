# palm_annot on RVMT reference sequences (per Phylum / Class / Order / Family).
# Run after 03_RDRP_identification.smk (through Step 22, split_by_phylum /
# split_by_rank) is complete -- {taxon}_ref.faa files must already exist in
# result/03_RDRP_identification/22_phylum_cluster/full_length_{rank}/.
#
# split_by_phylum / split_by_rank extract each taxon's top-hit RVMT reference
# sequences into {taxon}_ref.faa. Those references are already core-trimmed by
# RVMT's own pipeline, not by palmscan/palm_annot -- so to get a
# palmscan-confirmed A/B/C motif call (and eventually a palm_core/palm_extended
# region cut using the same boundary logic as extract_palm_regions.py) for the
# reference side, palm_annot needs to be run on them directly.
#
# Because this file is invoked as a separate Snakemake run after Step 22 has
# already finished, the {rank}/{taxon} wildcard values are read directly off
# disk at parse time below -- no checkpoint is needed (unlike a
# single-invocation pipeline where the producing rule and consuming rule are
# in the same DAG and the file doesn't exist yet when the DAG is first built).

import os

configfile: "config/config.yaml"

_PHYLUM_DIR      = "result/03_RDRP_identification/22_phylum_cluster"
_INPUT_DIR       = "result/03c_ref_palm_annot/00_input"
_REFPALM_DIR     = "result/03c_ref_palm_annot/1_palm_annot"
_RDRPREGION_DIR  = "result/03c_ref_palm_annot/2_rdrp_region"

# Which ranks to run palm_annot on -- set in config/config.yaml under
# ref_palm_annot.ranks (default: all 4). Restrict to e.g. ["phylum", "class"]
# to skip order/family for now.
_RANKS = config.get("ref_palm_annot", {}).get("ranks", ["phylum", "class", "order", "family"])

_REF_TAXA = {
    rank: sorted(
        f[:-len("_ref.faa")]
        for f in os.listdir(_PHYLUM_DIR + f"/full_length_{rank}")
        if f.endswith("_ref.faa")
    )
    for rank in _RANKS
}


rule all:
    input:
        [
            _REFPALM_DIR + f"/{rank}/{taxon}_ref_palmscan_hits.tsv"
            for rank in _RANKS
            for taxon in _REF_TAXA[rank]
        ],
        [
            _RDRPREGION_DIR + f"/{rank}/{taxon}/palm_regions.tsv"
            for rank in _RANKS
            for taxon in _REF_TAXA[rank]
        ],


# ── Step 00: stage a self-contained snapshot of the per-taxon ref fastas ─────

rule stage_ref_inputs:
    input:
        srcdir = _PHYLUM_DIR + "/full_length_{rank}",
    output:
        outdir = directory(_INPUT_DIR + "/full_length_{rank}"),
    wildcard_constraints:
        rank = "phylum|class|order|family",
    log:
        err = "log/03c_ref_palm_annot/00_input/stage_ref_inputs_{rank}.err",
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p {output.outdir} $(dirname {log.err})
        cp {input.srcdir}/*_ref.faa {output.outdir}/ 2> {log.err}
        """


rule palm_annot_ref:
    input:
        ref_fasta = _INPUT_DIR + "/full_length_{rank}/{taxon}_ref.faa",
    output:
        tsv     = _REFPALM_DIR + "/{rank}/{taxon}_ref_palmscan_hits.tsv",
        fev     = _REFPALM_DIR + "/{rank}/{taxon}_ref_palm_annot.fev",
        rdrp_aa = _REFPALM_DIR + "/{rank}/{taxon}_ref_rdrp_trimmed.faa",
    wildcard_constraints:
        rank = "phylum|class|order|family",
    log:
        log = "log/03c_ref_palm_annot/{rank}/{taxon}.log",
        err = "log/03c_ref_palm_annot/{rank}/{taxon}.err",
    params:
        palm_annot_py = config["palm_annot"]["install_dir"] + "/py/palm_annot.py",
        fev2tsv_py    = config["palm_annot"]["install_dir"] + "/py/fev2tsv.py",
        seqtype       = config["palm_annot"]["seqtype"],
        minscore      = config["palm_annot"]["minscore"],
        minpssmscore  = config["palm_annot"]["minpssmscore"],
        tmpdir        = config["palm_annot"]["tmpdir"] + "/ref_{rank}_{taxon}",
    threads: config["palm_annot"]["threads"]
    resources:
        mem_mb_per_cpu  = config["regular_memory"],
        runtime         = config["palm_annot"]["runtime"],
        cpus_per_task   = config["palm_annot"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account   = config["account"],
    conda:
        "../envs/palm_annot.yaml"
    shell:
        """
        mkdir -p $(dirname {output.tsv}) {params.tmpdir} $(dirname {log.log})
        python {params.palm_annot_py} \
            --input {input.ref_fasta} \
            --seqtype {params.seqtype} \
            --fev {output.fev} \
            --rdrp {output.rdrp_aa} \
            --minscore {params.minscore} \
            --minpssmscore {params.minpssmscore} \
            --threads {threads} \
            --tmpdir {params.tmpdir} \
            > {log.log} 2> {log.err}
        python {params.fev2tsv_py} \
            --input {output.fev} \
            --output {output.tsv} \
            --header yes \
            >> {log.log} 2>> {log.err}
        """


# ── Step 2: extract palm_core / palm_extended regions for each ref taxon ────
# Reuses the same script/logic as the acquired-sequence pipeline
# (scripts/03_RDRP_identification/extract_palm_regions.py): motif substrings
# (A/B/C) are read from rdrp_trimmed.faa's headers, then string-searched
# against the full (already core-trimmed by RVMT) reference sequence in
# {taxon}_ref.faa to compute core/extended boundaries. The --rdrpcatch-*
# flags are reused as-is (generic internally -- only affects an output
# subfolder name and a cosmetic "source" column value).

rule extract_ref_palm_regions:
    input:
        motif_fasta = _REFPALM_DIR + "/{rank}/{taxon}_ref_rdrp_trimmed.faa",
        full_fasta  = _INPUT_DIR + "/full_length_{rank}/{taxon}_ref.faa",
        script      = "scripts/03_RDRP_identification/extract_palm_regions.py",
    output:
        tsv      = _RDRPREGION_DIR + "/{rank}/{taxon}/palm_regions.tsv",
        core     = _RDRPREGION_DIR + "/{rank}/{taxon}/palm_core.faa",
        extended = _RDRPREGION_DIR + "/{rank}/{taxon}/palm_extended.faa",
        full     = _RDRPREGION_DIR + "/{rank}/{taxon}/rdrp_full.faa",
    wildcard_constraints:
        rank = "phylum|class|order|family",
    params:
        outdir     = _RDRPREGION_DIR + "/{rank}/{taxon}",
        rc_subdir  = _RDRPREGION_DIR + "/{rank}/{taxon}/RdRpCATCH",
        flank      = 150,
    log:
        out = "log/03c_ref_palm_annot/{rank}/{taxon}.extract.log",
        err = "log/03c_ref_palm_annot/{rank}/{taxon}.extract.err",
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
        mkdir -p {params.outdir} $(dirname {log.out})
        python {input.script} \
            --rdrpcatch-motif {input.motif_fasta} \
            --rdrpcatch-seqs  {input.full_fasta} \
            --outdir          {params.outdir} \
            --flank           {params.flank} \
            > {log.out} 2> {log.err}
        mv {params.rc_subdir}/* {params.outdir}/
        rmdir {params.rc_subdir}
        """
