import os

configfile: "config/config.yaml"

rule all:
    input:
        directory(config["bracken_report_dir"])


rule kraken2:
    input:
        r1 = config["rna_fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2 = config["rna_fastp_dir"] + "/{sample}/{sample}_2P.fq.gz"
    output:
        report = config["kraken2_dir"] + "/{sample}/{sample}_kraken2.report",
    log:
        log = config["kraken2_log_dir"] + "/kraken2/{sample}.log",
        err = config["kraken2_log_dir"] + "/kraken2/{sample}.err"
    conda:
        "../envs/kraken2.yaml"
    threads: config["kraken2"]["threads"]
    resources:
        mem_mb_per_cpu  = config["kraken2"]["memory"],
        runtime         = config["kraken2"]["runtime"],
        cpus_per_task   = config["kraken2"]["threads"],
        slurm_partition = config["kraken2"]["partition"],
        slurm_account   = config["kraken2"]["account"]
    params:
        out_dir    = config["kraken2_dir"] + "/{sample}",
        db         = config["kraken2_db"]["dir"],
        confidence = config["kraken2"]["confidence"]
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p $(dirname {log.log})
        kraken2 --db {params.db} \
                --threads {threads} \
                --confidence {params.confidence} \
                --report {output.report} \
                --paired {input.r1} {input.r2} \
                --output {output.out} \
                > {log.log} 2> {log.err}
        rm -f {output.out}
        """


rule bracken_build:
    input:
        db = directory(config["kraken2_db"]["dir"])
    output:
        marker = config["kraken2_db"]["dir"] + "/bracken_build.done"
    log:
        log = config["kraken2_log_dir"] + "/bracken_build/bracken_build.log",
        err = config["kraken2_log_dir"] + "/bracken_build/bracken_build.err"
    conda:
        "../envs/kraken2.yaml"
    threads: config["bracken_build"]["threads"]
    resources:
        mem_mb_per_cpu  = config["bracken_build"]["memory"],
        runtime         = config["bracken_build"]["runtime"],
        cpus_per_task   = config["bracken_build"]["threads"],
        slurm_partition = config["bracken_build"]["partition"],
        slurm_account   = config["bracken_build"]["account"]
    params:
        db          = config["kraken2_db"]["dir"],
        kmer_length = config["bracken_build"]["kmer_length"],
        read_length = config["bracken_build"]["read_length"]
    shell:
        """
        mkdir -p $(dirname {log.log})
        bracken-build -d {params.db} \
                      -t {threads} \
                      -k {params.kmer_length} \
                      -l {params.read_length} \
                      > {log.log} 2> {log.err}
        touch {output.marker}
        """


rule bracken:
    input:
        report     = config["kraken2_dir"] + "/{sample}/{sample}_kraken2.report",
        bracken_db = config["kraken2_db"]["dir"] + "/bracken_build.done"
    output:
        bracken_report_P = config["kraken2_dir"] + "/{sample}/{sample}_bracken_phylum.report",
        bracken_report_C = config["kraken2_dir"] + "/{sample}/{sample}_bracken_classes.report",
        bracken_report_O = config["kraken2_dir"] + "/{sample}/{sample}_bracken_orders.report",
        bracken_report_F = config["kraken2_dir"] + "/{sample}/{sample}_bracken_family.report",
        bracken_report_G = config["kraken2_dir"] + "/{sample}/{sample}_bracken_genus.report",
        bracken_report_S = config["kraken2_dir"] + "/{sample}/{sample}_bracken_species.report"
    log:
        log = config["kraken2_log_dir"] + "/bracken/{sample}.log",
        err = config["kraken2_log_dir"] + "/bracken/{sample}.err"
    conda:
        "../envs/kraken2.yaml"
    threads: config["bracken"]["threads"]
    resources:
        mem_mb_per_cpu  = config["bracken"]["memory"],
        runtime         = config["bracken"]["runtime"],
        cpus_per_task   = config["bracken"]["threads"],
        slurm_partition = config["bracken"]["partition"],
        slurm_account   = config["bracken"]["account"]
    params:
        db          = config["kraken2_db"]["dir"],
        read_length = config["bracken"]["read_length"]
    shell:
        """
        mkdir -p $(dirname {log.log})
        bracken -d {params.db} -i {input.report} -o {output.bracken_report_P} -r {params.read_length} -l P -t {threads} >> {log.log} 2>> {log.err}
        bracken -d {params.db} -i {input.report} -o {output.bracken_report_C} -r {params.read_length} -l C -t {threads} >> {log.log} 2>> {log.err}
        bracken -d {params.db} -i {input.report} -o {output.bracken_report_O} -r {params.read_length} -l O -t {threads} >> {log.log} 2>> {log.err}
        bracken -d {params.db} -i {input.report} -o {output.bracken_report_F} -r {params.read_length} -l F -t {threads} >> {log.log} 2>> {log.err}
        bracken -d {params.db} -i {input.report} -o {output.bracken_report_G} -r {params.read_length} -l G -t {threads} >> {log.log} 2>> {log.err}
        bracken -d {params.db} -i {input.report} -o {output.bracken_report_S} -r {params.read_length} -l S -t {threads} >> {log.log} 2>> {log.err}
        """


rule merge_bracken:
    input:
        P_reports = expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_phylum.report",  sample=config["rna_samples"]),
        C_reports = expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_classes.report",  sample=config["rna_samples"]),
        O_reports = expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_orders.report",   sample=config["rna_samples"]),
        F_reports = expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_family.report",   sample=config["rna_samples"]),
        G_reports = expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_genus.report",    sample=config["rna_samples"]),
        S_reports = expand(config["kraken2_dir"] + "/{sample}/{sample}_bracken_species.report",  sample=config["rna_samples"])
    output:
        directory(config["bracken_merge_dir"])
    threads: config["merge_bracken"]["threads"]
    resources:
        mem_mb_per_cpu  = config["merge_bracken"]["memory"],
        runtime         = config["merge_bracken"]["runtime"],
        cpus_per_task   = config["merge_bracken"]["threads"],
        slurm_partition = config["merge_bracken"]["partition"],
        slurm_account   = config["merge_bracken"]["account"]
    shell:
        """
        mkdir -p {output}/P {output}/C {output}/O {output}/F {output}/G {output}/S
        for f in {input.P_reports}; do cp $f {output}/P/; done
        for f in {input.C_reports}; do cp $f {output}/C/; done
        for f in {input.O_reports}; do cp $f {output}/O/; done
        for f in {input.F_reports}; do cp $f {output}/F/; done
        for f in {input.G_reports}; do cp $f {output}/G/; done
        for f in {input.S_reports}; do cp $f {output}/S/; done
        """


rule merge_bracken_report:
    input:
        dir = config["bracken_merge_dir"]
    output:
        dir = directory(config["bracken_report_dir"])
    conda:
        "../envs/python.yaml"
    threads: config["merge_bracken_report"]["threads"]
    resources:
        mem_mb_per_cpu  = config["merge_bracken_report"]["memory"],
        runtime         = config["merge_bracken_report"]["runtime"],
        cpus_per_task   = config["merge_bracken_report"]["threads"],
        slurm_partition = config["merge_bracken_report"]["partition"],
        slurm_account   = config["merge_bracken_report"]["account"]
    params:
        script = "scripts/101_kraken_check/merge_profiling_reports.py",
        outdir = lambda wildcards, output: os.path.abspath(output.dir)
    shell:
        """
        mkdir -p {params.outdir}
        python {params.script} -i {input.dir}/P -o {params.outdir}/phylum
        python {params.script} -i {input.dir}/C -o {params.outdir}/classes
        python {params.script} -i {input.dir}/O -o {params.outdir}/orders
        python {params.script} -i {input.dir}/F -o {params.outdir}/family
        python {params.script} -i {input.dir}/G -o {params.outdir}/genus
        python {params.script} -i {input.dir}/S -o {params.outdir}/species
        """

