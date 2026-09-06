#!/usr/bin/env python3
"""
Yangtze River RNA Virus Pipeline - Snakemake Workflow
RNA virus bioinformatics pipeline for processing RNA NGS sequencing data
"""

configfile: "config/config.yaml"


rule all:
    input:
        expand(
            config["rna_reformated_scaffolds_dir"] + "/rename_{min_len}/{sample}_scaffolds_rename_{min_len}.fasta",
            sample=config["rna_samples"],
            min_len=config["seqkit"]["min_length"]
        )


# Rule 1: Quality control and adapter removal with fastp
rule fastp_qc:
    input:
        r1 = ancient(config["rna_input_dir"] + "/{sample}/{sample}.R1.fq.gz"),
        r2 = ancient(config["rna_input_dir"] + "/{sample}/{sample}.R2.fq.gz")
    output:
        r1_paired   = config["rna_fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2_paired   = config["rna_fastp_dir"] + "/{sample}/{sample}_2P.fq.gz",
        r1_unpaired = config["rna_fastp_dir"] + "/{sample}/{sample}_U1.fq.gz",
        r2_unpaired = config["rna_fastp_dir"] + "/{sample}/{sample}_U2.fq.gz",
        html        = config["rna_fastp_dir"] + "/{sample}/{sample}.fastp.html",
        json        = config["rna_fastp_dir"] + "/{sample}/{sample}.fastp.json"
    log:
        log = config["rna_log_dir"] + "/fastp_qc/{sample}.log",
        err = config["rna_log_dir"] + "/fastp_qc/{sample}.err"
    params:
        quality_threshold = config["rna_fastp"]["quality_threshold"],
        length_required   = config["rna_fastp"]["length_required"]
    threads: config["rna_fastp"]["threads"]
    resources:
        mem_mb_per_cpu  = config["rna_fastp"]["memory_per_cpu"],
        runtime         = config["rna_fastp"]["runtime"],
        cpus_per_task   = config["rna_fastp"]["threads"],
        slurm_partition = config["rna_fastp"]["partition"],
        slurm_account   = config["rna_fastp"]["account"]
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        fastp --thread {threads} \
              --in1 {input.r1} --in2 {input.r2} \
              --out1 {output.r1_paired} --out2 {output.r2_paired} \
              --unpaired1 {output.r1_unpaired} --unpaired2 {output.r2_unpaired} \
              -h {output.html} -j {output.json} \
              --trim_poly_g --trim_poly_x \
              --qualified_quality_phred {params.quality_threshold} \
              --length_required {params.length_required} \
              > {log.log} 2> {log.err}
        """


# Rule 2: Remove rRNA sequences with ribodetector
rule ribodetector_rrna_removal:
    input:
        r1 = config["rna_fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2 = config["rna_fastp_dir"] + "/{sample}/{sample}_2P.fq.gz"
    output:
        r1_nonrrna = config["ribodetector_dir"] + "/{sample}/{sample}_nonrrna.1.fq.gz",
        r2_nonrrna = config["ribodetector_dir"] + "/{sample}/{sample}_nonrrna.2.fq.gz"
    log:
        log = config["rna_log_dir"] + "/ribodetector_rrna_removal/{sample}.log",
        err = config["rna_log_dir"] + "/ribodetector_rrna_removal/{sample}.err"
    params:
        min_length = config["ribodetector"]["min_length"],
        chunk_size = config["ribodetector"]["chunk_size"]
    threads: config["ribodetector"]["threads"]
    resources:
        mem_mb_per_cpu  = config["ribodetector"]["memory"],
        runtime         = config["ribodetector"]["runtime"],
        cpus_per_task   = config["ribodetector"]["threads"],
        slurm_partition = config["ribodetector"]["partition"],
        slurm_account   = config["ribodetector"]["account"]
    conda:
        "../envs/ribodetector.yaml"
    shell:
        """
        ribodetector_cpu -t {threads} \
                         -l {params.min_length} \
                         -i {input.r1} {input.r2} \
                         -e rrna \
                         --chunk_size {params.chunk_size} \
                         -o {output.r1_nonrrna} {output.r2_nonrrna} \
                         > {log.log} 2> {log.err}
        """


# Rule 3: Assembly with SPAdes
rule spades_assembly:
    input:
        r1 = config["ribodetector_dir"] + "/{sample}/{sample}_nonrrna.1.fq.gz",
        r2 = config["ribodetector_dir"] + "/{sample}/{sample}_nonrrna.2.fq.gz"
    output:
        scaffolds = config["rna_spades_dir"] + "/{sample}/{sample}_no_correction/scaffolds.fasta"
    log:
        log = config["rna_log_dir"] + "/spades_assembly/{sample}.log",
        err = config["rna_log_dir"] + "/spades_assembly/{sample}.err"
    params:
        outdir = config["rna_spades_dir"] + "/{sample}/{sample}_no_correction",
        k_list = config["rna_spades"]["k_values"],
        memory = lambda wildcards, resources: resources.mem_mb_per_cpu * resources.cpus_per_task // 1024
    threads: config["rna_spades"]["threads"]
    resources:
        mem_mb_per_cpu  = config["rna_spades"]["memory"],
        runtime         = config["rna_spades"]["runtime"],
        cpus_per_task   = config["rna_spades"]["threads"],
        slurm_partition = config["rna_spades"]["partition"],
        slurm_account   = config["rna_spades"]["account"]
    conda:
        "../envs/seqkit-spade.yaml"
    shell:
        """
        spades.py --meta \
            -o {params.outdir} \
            -1 {input.r1} -2 {input.r2} \
            -t {threads} -m {params.memory} \
            -k {params.k_list} \
            --only-assembler \
            > {log.log} 2> {log.err}
        """


# Rule 4: Rename and filter assemblies
rule rename_filter_assemblies:
    input:
        scaffolds = config["rna_spades_dir"] + "/{sample}/{sample}_no_correction/scaffolds.fasta"
    output:
        reformated = config["rna_reformated_scaffolds_dir"] + "/rename_{min_len}/{sample}_scaffolds_rename_{min_len}.fasta".format(min_len=config["seqkit"]["min_length"])
    log:
        log = config["rna_log_dir"] + "/rename_filter_assemblies/{sample}.log",
        err = config["rna_log_dir"] + "/rename_filter_assemblies/{sample}.err"
    params:
        min_length = config["seqkit"]["min_length"],
        nr_width   = config["seqkit"]["nr_width"]
    threads: config["seqkit"]["threads"]
    resources:
        mem_mb_per_cpu  = config["seqkit"]["memory"],
        runtime         = config["seqkit"]["runtime"],
        cpus_per_task   = config["seqkit"]["threads"],
        slurm_partition = config["seqkit"]["partition"],
        slurm_account   = config["seqkit"]["account"]
    conda:
        "../envs/seqkit-spade.yaml"
    shell:
        """
        seqkit seq -m {params.min_length} {input.scaffolds} 2> {log.err} | \
        seqkit replace -p .+ -r "{wildcards.sample}_{{nr}}" --nr-width {params.nr_width} \
        -o {output.reformated} \
        >> {log.log} 2>> {log.err}
        """
