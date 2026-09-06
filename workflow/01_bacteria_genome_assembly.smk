# this pipeline is used for bacteria genome assembly, From raw read , Fastp, and assembly by SPAdes, and rename contig as format structure: >Bac_{sample}_{10-digit-number}
configfile: "config/config.yaml"


rule all:
    input:
        expand(
            config["bacteria_reformated_scaffolds_dir"] + "/{sample}/{sample}_reformated.fa",
            sample=config["bacteria_samples"]
        ),
        # expand(
        #     config["bacteria_fastp_dir"] + "/{sample}/.cleaned",
        #     sample=config["bacteria_samples"]
        # )


rule fastp_quality_control:
    input:
        r1=ancient(config["bacteria_input_dir"] + "/{sample}.R1.raw.fastq.gz"),
        r2=ancient(config["bacteria_input_dir"] + "/{sample}.R2.raw.fastq.gz")
    output:
        r1_clean=config["bacteria_fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2_clean=config["bacteria_fastp_dir"] + "/{sample}/{sample}_2P.fq.gz",
        html=config["bacteria_fastp_dir"] + "/{sample}/{sample}.fastp.html",
        json=config["bacteria_fastp_dir"] + "/{sample}/{sample}.fastp.json"
    log:
        log=config["bacteria_log_dir"] + "/fastp_quality_control/{sample}.log",
        err=config["bacteria_log_dir"] + "/fastp_quality_control/{sample}.err"
    params:
        threads=config["bacteria_fastp"]["threads"],
        qualified_quality_phred=config["bacteria_fastp"]["qualified_quality_phred"],
        length_required=config["bacteria_fastp"]["length_required"],
        fastp_dir=config["bacteria_fastp_dir"]
    threads: config["bacteria_fastp"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["bacteria_fastp"]["runtime"],
        cpus_per_task=config["bacteria_fastp"]["threads"],
        slurm_partition=config["regular_partition"],
        slurm_account=config["account"]
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        mkdir -p {params.fastp_dir}/{wildcards.sample}
        fastp --thread {params.threads} \
              --in1 {input.r1} \
              --in2 {input.r2} \
              --out1 {output.r1_clean} \
              --out2 {output.r2_clean} \
              -h {output.html} \
              -j {output.json} \
              --trim_poly_g \
              --trim_poly_x \
              --qualified_quality_phred {params.qualified_quality_phred} \
              --length_required {params.length_required} \
              > {log.log} 2> {log.err}
        """
#        rm -f {input.r1} {input.r2}

rule spades_assembly_and_reformat:
    input:
        r1=config["bacteria_fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2=config["bacteria_fastp_dir"] + "/{sample}/{sample}_2P.fq.gz"
    output:
        reformated=config["bacteria_reformated_scaffolds_dir"] + "/{sample}/{sample}_reformated.fa",
        renaming=config["bacteria_reformated_scaffolds_dir"] + "/{sample}/{sample}_renaming.tsv",
        original=config["bacteria_reformated_scaffolds_dir"] + "/{sample}/{sample}_original_contigs.fasta"
    log:
        log=config["bacteria_log_dir"] + "/spades_assembly_and_reformat/{sample}.log",
        err=config["bacteria_log_dir"] + "/spades_assembly_and_reformat/{sample}.err"
    params:
        spades_threads=config["bacteria_spades"]["threads"],
        spades_memory=config["bacteria_spades"]["memory"],
        kmers=config["bacteria_spades"]["kmers"],
        spades_outdir=config["bacteria_spades_dir"] + "/{sample}",
        spades_dir=config["bacteria_spades_dir"],
        min_length=config["reformat_scaffolds"]["min_length"],
        nr_width=config["reformat_scaffolds"]["nr_width"],
        reformated_scaffolds_dir=config["bacteria_reformated_scaffolds_dir"]
    threads: config["bacteria_spades"]["threads"]
    resources:
        mem_mb_per_cpu=config["bacteria_spades"]["memory_per_cpu"],
        runtime=config["bacteria_spades"]["runtime"] + config["reformat_scaffolds"]["runtime"],
        cpus_per_task=config["bacteria_spades"]["threads"],
        slurm_partition=config["bacteria_spades"]["partition"],
        slurm_account=config["account"]
    conda:
        "../envs/seqkit-spade.yaml"
    shell:
        """
        mkdir -p {params.spades_dir}/{wildcards.sample}
        spades.py --meta \
                  -o {params.spades_outdir} \
                  -1 {input.r1} \
                  -2 {input.r2} \
                  -t {params.spades_threads} \
                  -m {params.spades_memory} \
                  -k {params.kmers} \
                  --only-assembler > {log.log} 2> {log.err}
        mkdir -p {params.reformated_scaffolds_dir}/{wildcards.sample}
        cp {params.spades_outdir}/contigs.fasta {output.original}
        seqkit seq -m {params.min_length} {output.original} | \
        seqkit replace -p .+ -r "Bac_{wildcards.sample}_{{nr}}" --nr-width {params.nr_width} \
        -o {output.reformated} >> {log.log} 2>> {log.err}
        paste <(seqkit seq -m {params.min_length} {output.original} | seqkit fx2tab -n) \
              <(seqkit fx2tab -n {output.reformated}) \
        > {output.renaming}
        rm -rf {params.spades_dir}/{wildcards.sample}
        """

# rule cleanup_fastp:
#     input:
#         config["bacteria_reformated_scaffolds_dir"] + "/{sample}/{sample}_reformated.fa",
#     output:
#         touch(config["bacteria_fastp_dir"] + "/{sample}/.cleaned"),
#     log:
#         err=config["bacteria_log_dir"] + "/cleanup_fastp/{sample}.err",
#     threads: 1
#     resources:
#         mem_mb_per_cpu  = config["small_job"]["memory"],
#         runtime         = 30,
#         cpus_per_task   = 1,
#         slurm_partition = config["regular_partition"],
#         slurm_account   = config["account"],
#     shell:
#         """
#         rm -f {config[bacteria_fastp_dir]}/{wildcards.sample}/{wildcards.sample}_1P.fq.gz \
#               {config[bacteria_fastp_dir]}/{wildcards.sample}/{wildcards.sample}_2P.fq.gz \
#               2> {log.err}
#         """


# Step 3: Filter contigs >100bp and rename per sample
rule seqkit_overXbp:
    input:
        config["bacteria_reformated_scaffolds_dir"] + "/{sample}/{sample}_original_contigs.fasta"
    output:
        config["bacteria_reformated_scaffolds_dir"] + "/{sample}/{sample}_over" + str(config["seqkit_overXbp"]["min_length"]) + "bp.fa"
    log:
        log=config["bacteria_log_dir"] + "/seqkit_overXbp/{sample}.log",
        err=config["bacteria_log_dir"] + "/seqkit_overXbp/{sample}.err"
    params:
        min_length=config["seqkit_overXbp"]["min_length"],
        nr_width=config["seqkit_overXbp"]["nr_width"]
    threads: config["seqkit_overXbp"]["threads"]
    resources:
        mem_mb_per_cpu=config["regular_memory"],
        runtime=config["seqkit_overXbp"]["runtime"],
        cpus_per_task=config["seqkit_overXbp"]["threads"],
        slurm_partition=config["seqkit_overXbp"]["partition"],
        slurm_account=config["account"]
    conda:
        "../envs/seqkit-spade.yaml"
    shell:
        """
        seqkit seq -m {params.min_length} {input} | \
        seqkit replace -p .+ -r "Bac_{wildcards.sample}_{{nr}}" --nr-width {params.nr_width} \
        -o {output} > {log.log} 2> {log.err}
        """


# Step 4: Concatenate all samples into one file
rule cat_all_scaffolds:
    input:
        expand(
            config["bacteria_reformated_scaffolds_dir"] + "/{sample}/{sample}_over" + str(config["seqkit_overXbp"]["min_length"]) + "bp.fa",
            sample=config["bacteria_samples"]
        )
    output:
        config["easy-linclust_dir"] + "/all_samples_over100bp.fa"
    log:
        err=config["bacteria_log_dir"] + "/cat_all_scaffolds/cat_all_scaffolds.err"
    shell:
        "cat {input} > {output} 2> {log.err}"


# Step 5: Cluster all scaffolds with mmseqs easy-linclust to remove redundancy
rule mmseqs_easy_linclust:
    input:
        config["easy-linclust_dir"] + "/all_samples_over100bp.fa"
    output:
        rep_seq=config["easy-linclust_dir"] + "/linclust_rep_seq.fasta"
    log:
        log=config["bacteria_log_dir"] + "/bac_mmseqs_easy_linclust/mmseqs_easy_linclust.log",
        err=config["bacteria_log_dir"] + "/bac_mmseqs_easy_linclust/mmseqs_easy_linclust.err"
    params:
        prefix=config["easy-linclust_dir"] + "/linclust",
        tmp_dir=config["easy-linclust_dir"] + "/tmp",
        min_seq_id=config["easy-linclust"]["min_seq_id"],
        cov_mode=config["easy-linclust"]["cov_mode"],
        cluster_mode=config["easy-linclust"]["cluster_mode"],
        coverage=config["easy-linclust"]["coverage"]
    threads: config["easy-linclust"]["threads"]
    resources:
        mem_mb_per_cpu=config["easy-linclust"]["memory_per_cpu"],
        runtime=config["easy-linclust"]["runtime"],
        cpus_per_task=config["easy-linclust"]["threads"],
        slurm_partition=config["easy-linclust"]["partition"],
        slurm_account=config["account"]
    conda:
        "../envs/mmseqs.yaml"
    shell:
        """
        mmseqs easy-linclust {input} {params.prefix} {params.tmp_dir} \
            --min-seq-id {params.min_seq_id} \
            --cov-mode {params.cov_mode} \
            --cluster-mode {params.cluster_mode} \
            -c {params.coverage} \
            --threads {threads} \
            > {log.log} 2> {log.err}
        """
