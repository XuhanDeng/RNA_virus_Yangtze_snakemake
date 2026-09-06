configfile: "config/config.yaml"

SAMPLES    = config["rna_samples"]
_LOG_DIR   = "log/08_other_virus_identification"
_RESULT    = "result/08_other_virus_identification"
_FASTA_DIR = config["rna_reformated_scaffolds_dir"]
_MIN_LEN   = str(config["seqkit"]["min_length"])
_FASTA     = lambda wc: (
    _FASTA_DIR + f"/rename_{_MIN_LEN}/{wc.sample}_scaffolds_rename_{_MIN_LEN}.fasta"
)

_MERGED_SUMMARY    = _RESULT + "/2_virus_filter_summary/all_samples_virus_filter_summary.tsv"
_ANNOTATED_SUMMARY = _RESULT + "/2_virus_filter_summary/all_samples_virus_filter_summary.annotated.tsv"


rule all:
    input:
        expand(
            _RESULT + "/1_genomad/{sample}/{sample}_scaffolds_rename_" + _MIN_LEN + "_summary/{sample}_scaffolds_rename_" + _MIN_LEN + "_virus.fna",
            sample=SAMPLES
        ),
        expand(_RESULT + "/1_genomad/{sample}/{sample}_virus_filter_summary.tsv", sample=SAMPLES),
        _ANNOTATED_SUMMARY


# ── Step 1: GeNomad identification (per sample) ──────────────────────────────

rule genomad_identification:
    input:
        fasta = _FASTA
    output:
        fna = _RESULT + "/1_genomad/{sample}/{sample}_scaffolds_rename_" + _MIN_LEN + "_summary/{sample}_scaffolds_rename_" + _MIN_LEN + "_virus.fna"
    conda:
        "../envs/genomad.yaml"
    log:
        log = _LOG_DIR + "/1_genomad/{sample}.log",
        err = _LOG_DIR + "/1_genomad/{sample}.err"
    threads: config["genomad"]["threads"]
    resources:
        slurm_partition = config["genomad"]["partition"],
        runtime         = config["genomad"]["runtime"],
        mem_mb_per_cpu  = config["genomad"]["memory"],
        cpus_per_task   = config["genomad"]["threads"],
        slurm_account   = config["genomad"]["account"]
    params:
        db     = config["genomad"]["db_dir"],
        outdir = _RESULT + "/1_genomad/{sample}"
    shell:
        """
        mkdir -p {params.outdir}
        mkdir -p {_LOG_DIR}/1_genomad
        genomad end-to-end --cleanup -t {threads} \
            {input.fasta} \
            {params.outdir} \
            {params.db} \
            > {log.log} 2> {log.err}
        """


# ── Step 2: Filter + merge + annotate ────────────────────────────────────────

rule filter_virus_summary:
    input:
        summary = _RESULT + "/1_genomad/{sample}/{sample}_scaffolds_rename_" + _MIN_LEN + "_summary/{sample}_scaffolds_rename_" + _MIN_LEN + "_virus_summary.tsv"
    output:
        tsv = _RESULT + "/1_genomad/{sample}/{sample}_virus_filter_summary.tsv"
    conda:
        "../envs/python.yaml"
    log:
        err = _LOG_DIR + "/1_genomad/{sample}_virus_filter.err"
    threads: config["small_job"]["threads"]
    resources:
        slurm_partition = config["small_job"]["partition"],
        runtime         = config["small_job"]["runtime"],
        mem_mb_per_cpu  = config["small_job"]["memory"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_account   = config["small_job"]["account"]
    params:
        exclude_topology = config["genomad_filter"]["exclude_topology"],
        min_hallmarks    = config["genomad_filter"]["min_hallmarks"]
    shell:
        """
        python scripts/08_other_virus_identification/filter_virus_summary.py \
            --summary          {input.summary} \
            --output           {output.tsv} \
            --exclude-topology {params.exclude_topology} \
            --min-hallmarks    {params.min_hallmarks} \
            2> {log.err}
        """


# ── Step 3: Merge all-sample filtered summaries ───────────────────────────────

rule merge_virus_summary:
    input:
        expand(_RESULT + "/1_genomad/{sample}/{sample}_virus_filter_summary.tsv", sample=SAMPLES)
    output:
        _MERGED_SUMMARY
    log:
        err = _LOG_DIR + "/2_virus_filter_summary.err"
    threads: config["small_job"]["threads"]
    resources:
        slurm_partition = config["small_job"]["partition"],
        runtime         = config["small_job"]["runtime"],
        mem_mb_per_cpu  = config["small_job"]["memory"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_account   = config["small_job"]["account"]
    shell:
        """
        mkdir -p $(dirname {output})
        head -1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f"; done >> {output} 2> {log.err}
        """


# ── Step 4: Annotate with ICTV VMR taxonomy + genome composition ──────────────

rule annotate_with_ICTV_VMR:
    input:
        summary = _MERGED_SUMMARY,
        vmr     = config["ICTV_VMR"]["xlsx"]
    output:
        _ANNOTATED_SUMMARY
    conda:
        "../envs/python.yaml"
    log:
        err = _LOG_DIR + "/2_annotate_ICTV_VMR.err"
    threads: config["small_job"]["threads"]
    resources:
        slurm_partition = config["small_job"]["partition"],
        runtime         = config["small_job"]["runtime"],
        mem_mb_per_cpu  = config["small_job"]["memory"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_account   = config["small_job"]["account"]
    params:
        sheet = config["ICTV_VMR"]["sheet"]
    shell:
        """
        python scripts/08_other_virus_identification/annotate_ICTV_VMR.py \
            --summary {input.summary} \
            --vmr     {input.vmr} \
            --sheet   "{params.sheet}" \
            --output  {output} \
            2> {log.err}
        """
