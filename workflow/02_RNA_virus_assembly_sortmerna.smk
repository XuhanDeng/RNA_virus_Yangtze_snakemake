#!/usr/bin/env python3
"""
Yangtze River RNA Virus Pipeline - Snakemake Workflow
RNA virus bioinformatics pipeline for processing RNA NGS sequencing data
"""

configfile: "config/config.yaml"


rule all:
    input:
        expand(
            config["rna_reformated_scaffolds_dir"] + "/rename_1000/{sample}_scaffolds_rename_1000.fasta",
            sample=config["rna_samples"]
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
        slurm_account   = config["account"]
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


# Rule 2: Remove rRNA sequences with SortMeRNA
rule sortmerna_rrna_removal:
    input:
        r1 = config["rna_fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2 = config["rna_fastp_dir"] + "/{sample}/{sample}_2P.fq.gz"
    output:
        r1_nonrrna = config["sortmerna_dir"] + "/{sample}/{sample}_nonrrna_fwd.fq.gz",
        r2_nonrrna = config["sortmerna_dir"] + "/{sample}/{sample}_nonrrna_rev.fq.gz"
    log:
        log = config["rna_log_dir"] + "/sortmerna_rrna_removal/{sample}.log",
        err = config["rna_log_dir"] + "/sortmerna_rrna_removal/{sample}.err"
    params:
        other_prefix = config["sortmerna_dir"] + "/{sample}/{sample}_nonrrna",
        workdir      = config["sortmerna_dir"] + "/{sample}/workdir",
        ref_flags    = lambda wildcards: " ".join(
                           "--ref " + db for db in config["sortmerna"]["dbs"]
                       )
    threads: config["sortmerna"]["threads"]
    resources:
        mem_mb_per_cpu  = config["sortmerna"]["memory"],
        runtime         = config["sortmerna"]["runtime"],
        cpus_per_task   = config["sortmerna"]["threads"],
        slurm_partition = config["sortmerna"]["partition"],
        slurm_account   = config["sortmerna"]["account"]
    conda:
        "../envs/sortmerna.yaml"
    shell:
        """
        mkdir -p {params.workdir}
        sortmerna \
            {params.ref_flags} \
            --reads {input.r1} --reads {input.r2} \
            --other {params.other_prefix} \
            --fastx \
            --out2 \
            --paired-in \
            --threads 10 \
            --workdir {params.workdir} \
            > {log.log} 2> {log.err}
        """
'''
remeber change --threas 10 to {threads} after testing
'''

# Rule 3: Re-pair reads after SortMeRNA with BBTools repair.sh
rule bbtools_repair:
    input:
        r1 = config["sortmerna_dir"] + "/{sample}/{sample}_nonrrna_fwd.fq.gz",
        r2 = config["sortmerna_dir"] + "/{sample}/{sample}_nonrrna_rev.fq.gz"
    output:
        r1 = config["bbtools_dir"] + "/{sample}/{sample}_paired_1.fq.gz",
        r2 = config["bbtools_dir"] + "/{sample}/{sample}_paired_2.fq.gz"
    log:
        log = config["rna_log_dir"] + "/bbtools_repair/{sample}.log",
        err = config["rna_log_dir"] + "/bbtools_repair/{sample}.err"
    conda:
        "../envs/bbtools.yaml"
    threads: config["bbtools_repair"]["threads"]
    resources:
        mem_mb_per_cpu  = config["bbtools_repair"]["memory"],
        runtime         = config["bbtools_repair"]["runtime"],
        cpus_per_task   = config["bbtools_repair"]["threads"],
        slurm_partition = config["bbtools_repair"]["partition"],
        slurm_account   = config["bbtools_repair"]["account"]
    params:
        out_dir = config["bbtools_dir"] + "/{sample}"
    shell:
        """
        mkdir -p {params.out_dir}
        repair.sh \
            in1={input.r1} in2={input.r2} \
            out1={output.r1} out2={output.r2} \
            repair=t \
            > {log.log} 2> {log.err}
        """


# Rule 4: Assembly with SPAdes
rule spades_assembly:
    input:
        r1 = config["bbtools_dir"] + "/{sample}/{sample}_paired_1.fq.gz",
        r2 = config["bbtools_dir"] + "/{sample}/{sample}_paired_2.fq.gz"
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
        reformated = config["rna_reformated_scaffolds_dir"] + "/rename_1000/{sample}_scaffolds_rename_1000.fasta"
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
