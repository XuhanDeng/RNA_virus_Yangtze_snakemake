#!/usr/bin/env python3
"""
Yangtze River RNA Virus Pipeline - Snakemake Workflow
RNA virus bioinformatics pipeline for processing RNA NGS sequencing data
"""

# Configuration
configfile: "config/config.yaml"


# Helper to build optional CLI args for Ref_cluster_select.py
def build_ref_cluster_args(cfg):
    params = cfg.get("ref_cluster_select", {})
    args = []
    db_id_col = params.get("db_id_col")
    if db_id_col:
        args += ["--db-id-col", str(db_id_col)]
    sample_order = params.get("sample_order")
    if sample_order:
        args += ["--sample-order", str(sample_order)]
    sample_order_file = params.get("sample_order_file")
    if sample_order_file:
        args += ["--sample-order-file", str(sample_order_file)]
    if params.get("with_cluster_id"):
        args.append("--with-cluster-id")
    return " ".join(args)

REF_CLUSTER_ARGS = build_ref_cluster_args(config)


# Final output rule
rule all:
    input:
        "result/04_modified_esvirtue/databases/final_merged_database_only_assembly/virus_pathogen_database.all_metadata.tsv",
        "result/04_modified_esvirtue/databases/final_merged_database_only_assembly/virus_pathogen_database.fna",
        "result/04_modified_esvirtue/databases/final_merged_database_only_assembly/virus_pathogen_database.mmi",
        expand("result/04_modified_esvirtue/es/first_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv", sample=config["rna_samples"]),
        expand("result/04_modified_esvirtue/es/second_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv", sample=config["rna_samples"]),
        expand("result/04_modified_esvirtue/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.rpkmf.tsv", sample=config["rna_samples"]),
        "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.tsv",
        "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.subspecies.tsv",
        "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.read_count.tsv",
        "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.subspecies.read_count.tsv",
        "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.species.read_count.tsv",





# Step5.1 Merge with Esviritu database_just_merge
# Esviritu is a tool for identifying viral sequences from metagenomic data. however, After I try to mapping my filtered sequencing data to Esvirtu database, I found very few reads can be mapped to the database. So I decide to enlarge the database by merging established Esviritu database with my assembled contigs from metaspades. The minimum contig length is set to 200bp to ensure quality. which aligns with the Esviritu database construction criteria. for this step i just use seqkit to filter and cat all database together.

rule merge_esviritu_database_1:
    input:
        contigs = expand(
            config["rna_reformated_scaffolds_dir"] + "/rename_{min_len}/{sample}_scaffolds_rename_{min_len}.fasta",
            sample=config["rna_samples"],
            min_len=config["seqkit"]["min_length"],
        )
    output:
        merged_db = "result/04_modified_esvirtue/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta"
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["merge_esviritu_database_1"]["threads"]
    resources:
        mem_mb_per_cpu  = config["merge_esviritu_database_1"]["memory"],
        runtime         = config["merge_esviritu_database_1"]["runtime"],
        cpus_per_task   = config["merge_esviritu_database_1"]["threads"],
        slurm_partition = config["merge_esviritu_database_1"]["partition"],
        slurm_account   = config["merge_esviritu_database_1"]["account"]
    params:
        min_length = config["merge_esviritu_database_1"]["min_length"]
    log:
        out="log/04_modified_esvirtue/merge_esviritu_database/merge_esviritu_database_only_assembly_1.log",
        err="log/04_modified_esvirtue/merge_esviritu_database/merge_esviritu_database_only_assembly_1.err"
    shell:
        """
        mkdir -p result/04_modified_esvirtue/databases/database_merged_only_assembly
        seqkit seq -m {params.min_length} {input.contigs} >> {output.merged_db}
        """

# Previous CheckV-based clustering strategy (kept for reference)
rule merge_esviritu_database_2_checkv:
    input:
        merged_db="result/04_modified_esvirtue/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta"
    output:
        blast_db=directory("result/04_modified_esvirtue/databases/cluster_only_assembly/blast_db"),
        blast_results="result/04_modified_esvirtue/databases/cluster_only_assembly/all_samples_blast.tsv",
        ani_results="result/04_modified_esvirtue/databases/cluster_only_assembly/all_samples_ani.tsv",
        cluster_results="result/04_modified_esvirtue/databases/cluster_only_assembly/all_samples_cluster.tsv"
    conda: "../envs/checkv.yaml"
    threads: config["checkv"]["threads"]
    params:
        blast_outfmt        = config["checkv"]["blast_outfmt"],
        blast_max_target    = config["checkv"]["blast_max_target_seqs"],
        min_ani             = config["checkv"]["min_ani"],
        min_coverage        = config["checkv"]["min_coverage"],
        min_qcov            = config["checkv"]["min_qcov"],
        checkv_ani          = config["scripts"]["checkv_ani"],
        checkv_clust        = config["scripts"]["checkv_clust"]
    log:
        out="log/04_modified_esvirtue/cluster_all.log",
        err="log/04_modified_esvirtue/cluster_all.err"
    resources:
        mem_mb_per_cpu  = config["checkv"]["memory"],
        runtime         = config["checkv"]["runtime"],
        cpus_per_task   = config["checkv"]["threads"],
        slurm_partition = config["checkv"]["partition"],
        slurm_account   = config["checkv"]["account"]
    shell:
        """
        mkdir -p result/04_modified_esvirtue/databases/cluster_only_assembly
        mkdir -p log/04_modified_esvirtue

        # Create BLAST database
        makeblastdb -in {input.merged_db} \
            -dbtype nucl \
            -out result/04_modified_esvirtue/databases/cluster_only_assembly/all_samples_db \
            >> {log.out} 2>> {log.err}

        # Run BLAST all-vs-all
        blastn -query {input.merged_db} \
            -db result/04_modified_esvirtue/databases/cluster_only_assembly/all_samples_db \
            -outfmt '{params.blast_outfmt}' \
            -max_target_seqs {params.blast_max_target} \
            -out {output.blast_results} \
            -num_threads {threads} \
            >> {log.out} 2>> {log.err}

        # Calculate ANI
        python {params.checkv_ani} \
            -i {output.blast_results} \
            -o {output.ani_results} \
            >> {log.out} 2>> {log.err}

        # Cluster sequences
        python {params.checkv_clust} \
            --fna {input.merged_db} \
            --ani {output.ani_results} \
            --out {output.cluster_results} \
            --min_ani {params.min_ani} \
            --min_tcov {params.min_coverage} \
            --min_qcov {params.min_qcov} \
            >> {log.out} 2>> {log.err}

        # Create directory for blast database files
        mkdir -p {output.blast_db}
        mv result/04_modified_esvirtue/databases/cluster_only_assembly/all_samples_db.* {output.blast_db}/
        """


# Select representative sequences from each cluster
rule select_cluster_representatives:
    input:
        cluster="result/04_modified_esvirtue/databases/cluster_only_assembly/all_samples_cluster.tsv",
    output:
        reps="result/04_modified_esvirtue/databases/cluster_only_assembly/cluster_representatives.txt",
    conda:
        "../envs/python.yaml"
    threads: config["select_cluster_representatives"]["threads"]
    resources:
        mem_mb_per_cpu  = config["select_cluster_representatives"]["memory"],
        runtime         = config["select_cluster_representatives"]["runtime"],
        cpus_per_task   = config["select_cluster_representatives"]["threads"],
        slurm_partition = config["select_cluster_representatives"]["partition"],
        slurm_account   = config["select_cluster_representatives"]["account"]
    params:
        script = config["scripts"]["ref_cluster_select_assembly_only"]
    log:
        out="log/04_modified_esvirtue/cluster_representatives.log",
        err="log/04_modified_esvirtue/cluster_representatives.err",
    shell:
        """
        mkdir -p log/04_modified_esvirtue
        python {params.script} \
            --cluster {input.cluster} \
            --output {output.reps} \
            {REF_CLUSTER_ARGS} \
            > {log.out} 2> {log.err}
        """

rule extract_final_virus_pathogen_database:
    input:
        reps="result/04_modified_esvirtue/databases/cluster_only_assembly/cluster_representatives.txt",
        fasta="result/04_modified_esvirtue/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta",
    output:
        fasta="result/04_modified_esvirtue/databases/final_merged_database_only_assembly/virus_pathogen_database.fna",
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        mem_mb_per_cpu  = config["seqkit"]["memory"],
        runtime         = config["seqkit"]["runtime"],
        cpus_per_task   = config["seqkit"]["threads"],
        slurm_partition = config["seqkit"]["partition"],
        slurm_account   = config["seqkit"]["account"]
    log:
        out="log/04_modified_esvirtue/extract_final_db.log",
        err="log/04_modified_esvirtue/extract_final_db.err"
    shell:
        """
        mkdir -p result/04_modified_esvirtue/databases/final_merged_database_only_assembly
        mkdir -p log/04_modified_esvirtue
        seqkit grep -f {input.reps} {input.fasta} -o {output.fasta} \
            > {log.out} 2> {log.err}
        """

rule index_final_virus_pathogen_database:
    input:
        fasta="result/04_modified_esvirtue/databases/final_merged_database_only_assembly/virus_pathogen_database.fna",
    output:
        mmi="result/04_modified_esvirtue/databases/final_merged_database_only_assembly/virus_pathogen_database.mmi",
    conda:
        "../envs/minimap2.yaml"
    threads: config["index_final_virus_pathogen_database"]["threads"]
    resources:
        mem_mb_per_cpu  = config["index_final_virus_pathogen_database"]["memory"],
        runtime         = config["index_final_virus_pathogen_database"]["runtime"],
        cpus_per_task   = config["index_final_virus_pathogen_database"]["threads"],
        slurm_partition = config["index_final_virus_pathogen_database"]["partition"],
        slurm_account   = config["index_final_virus_pathogen_database"]["account"]
    log:
        out="log/04_modified_esvirtue/minimap2_index.log",
        err="log/04_modified_esvirtue/minimap2_index.err"
    shell:
        """
        mkdir -p result/04_modified_esvirtue/databases/final_merged_database_only_assembly
        mkdir -p log/04_modified_esvirtue
        minimap2 -d {output.mmi} {input.fasta} \
            > {log.out} 2> {log.err}
        """

rule seqkit_length_from_merged_db:
    input:
        fasta="result/04_modified_esvirtue/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta",
    output:
        length="result/04_modified_esvirtue/databases/final_merged_database_only_assembly/length.txt",
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        mem_mb_per_cpu  = config["seqkit"]["memory"],
        runtime         = config["seqkit"]["runtime"],
        cpus_per_task   = config["seqkit"]["threads"],
        slurm_partition = config["seqkit"]["partition"],
        slurm_account   = config["seqkit"]["account"]
    log:
        out="log/04_modified_esvirtue/seqkit_length.log",
        err="log/04_modified_esvirtue/seqkit_length.err"
    shell:
        """
        mkdir -p result/04_modified_esvirtue/databases/final_merged_database_only_assembly
        mkdir -p log/04_modified_esvirtue
        seqkit fx2tab -n -l {input.fasta} > {output.length} \
            2> {log.err}
        """

rule merge_esviritu_metadata:
    input:
        id_length="result/04_modified_esvirtue/databases/final_merged_database_only_assembly/length.txt",
        ids="result/04_modified_esvirtue/databases/cluster_only_assembly/cluster_representatives.txt",
        fasta="result/04_modified_esvirtue/databases/database_merged_only_assembly/esviritu_merged_db_only_assembly.fasta",
    output:
        merged="result/04_modified_esvirtue/databases/final_merged_database_only_assembly/virus_pathogen_database.all_metadata.tsv",
    conda:
        "../envs/python.yaml"
    threads: config["merge_esviritu_metadata"]["threads"]
    resources:
        mem_mb_per_cpu  = config["merge_esviritu_metadata"]["memory"],
        runtime         = config["merge_esviritu_metadata"]["runtime"],
        cpus_per_task   = config["merge_esviritu_metadata"]["threads"],
        slurm_partition = config["merge_esviritu_metadata"]["partition"],
        slurm_account   = config["merge_esviritu_metadata"]["account"]
    params:
        script = config["scripts"]["merge_es_tsv_assembly_only"]
    log:
        out="log/04_modified_esvirtue/merge_metadata.log",
        err="log/04_modified_esvirtue/merge_metadata.err"
    shell:
        """
        mkdir -p result/04_modified_esvirtue/databases/final_merged_database_only_assembly
        mkdir -p log/04_modified_esvirtue
        python {params.script} \
            --id-length {input.id_length} \
            --ids {input.ids} \
            --fasta {input.fasta} \
            --output {output.merged} \
            > {log.out} 2> {log.err}
        """

































rule esviritu_mapped_to_standard_database:
    input:
        r1=config["rna_fastp_dir"] + "/{sample}/{sample}_1P.fq.gz",
        r2=config["rna_fastp_dir"] + "/{sample}/{sample}_2P.fq.gz",
        db_dir=config["databases"]["esviritu_db_v3.2.4"],
    output:
        assembly_summary="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv",
        info="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}.detected_virus.info.tsv",
        consensus="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_final_consensus.fasta",
        coverage="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}.virus_coverage_windows.tsv",
        log_file="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_esviritu.log",
        params="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_esviritu.params.yaml",
        tagged_r1="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_nonrrna.1.tagged.fq.gz",
        tagged_r2="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_nonrrna.2.tagged.fq.gz",
        mapped_reads="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_temp/{sample}.reads.txt",
    conda:
        "../envs/esviritu_map.yaml"
    threads: config["esviritu_map"]["threads"]
    resources:
        mem_mb_per_cpu  = config["esviritu_map"]["memory"],
        runtime         = config["esviritu_map"]["runtime"],
        cpus_per_task   = config["esviritu_map"]["threads"],
        slurm_partition = config["esviritu_map"]["partition"],
        slurm_account   = config["esviritu_map"]["account"]
    params:
        temp_dir    = "result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_temp",
        third_bam   = "result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_temp/{sample}.third.filt.sorted.bam",
        keep        = config["esviritu_map"]["keep"],
        quiet       = config["esviritu_map"]["quiet"]
    log:
        out="log/04_modified_esvirtue/first_filter/{sample}.log",
        err="log/04_modified_esvirtue/first_filter/{sample}.err"
    shell:
        """
        mkdir -p result/04_modified_esvirtue/es/first_filter/{wildcards.sample}
        mkdir -p log/04_modified_esvirtue/first_filter
        zcat {input.r1} | gawk 'NR%4==1{{ $0=gensub(/^@([^ ]+)/, "@\\\\1_1", 1) }} {{ print }}' | gzip > {output.tagged_r1}
        zcat {input.r2} | gawk 'NR%4==1{{ $0=gensub(/^@([^ ]+)/, "@\\\\1_2", 1) }} {{ print }}' | gzip > {output.tagged_r2}
        EsViritu -r {output.tagged_r1} {output.tagged_r2} \
                 -s {wildcards.sample} \
                 -t {threads} \
                 -o result/04_modified_esvirtue/es/first_filter/{wildcards.sample} \
                 --db {input.db_dir} \
                 --keep {params.keep} -q {params.quiet} \
                 > {log.out} 2>> {log.err}
        mkdir -p {params.temp_dir}
        samtools view {params.third_bam} | cut -f1 | sort -u > {output.mapped_reads}
        """


rule extract_unmapped_reads_esviritu:
    input:
        r1="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_nonrrna.1.tagged.fq.gz",
        r2="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_nonrrna.2.tagged.fq.gz",
        mapped_reads="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_temp/{sample}.reads.txt",
    output:
        unmapped_r1="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_nonrrna.1.unmapped.fq.gz",
        unmapped_r2="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_nonrrna.2.unmapped.fq.gz",
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        mem_mb_per_cpu  = config["seqkit"]["memory"],
        runtime         = config["seqkit"]["runtime"],
        cpus_per_task   = config["seqkit"]["threads"],
        slurm_partition = config["seqkit"]["partition"],
        slurm_account   = config["seqkit"]["account"]
    log:
        out="log/04_modified_esvirtue/unmapped_reads/{sample}.log",
        err="log/04_modified_esvirtue/unmapped_reads/{sample}.err"
    shell:
        """
        mkdir -p result/04_modified_esvirtue/es/first_filter/{wildcards.sample}
        mkdir -p log/04_modified_esvirtue/unmapped_reads
        tmp_dir=result/04_modified_esvirtue/es/first_filter/{wildcards.sample}/{wildcards.sample}_temp
        mkdir -p $tmp_dir
        seqkit seq -n {input.r1} > $tmp_dir/fastq_full_ids.txt
        seqkit seq -n {input.r2} >> $tmp_dir/fastq_full_ids.txt
        awk '
          NR==FNR {{
            core=$1
            full=$0
            map[core]=full
            next
          }}
          {{
            core=$1
            if (core in map) print map[core]
            else print core "\\tNOT_FOUND" > "/dev/stderr"
          }}
        ' $tmp_dir/fastq_full_ids.txt {input.mapped_reads} > $tmp_dir/ids.full.txt 2> {log.err}
        seqkit grep -v -n -f $tmp_dir/ids.full.txt {input.r1} -o {output.unmapped_r1} \
            >> {log.out} 2>> {log.err}
        seqkit grep -v -n -f $tmp_dir/ids.full.txt {input.r2} -o {output.unmapped_r2} \
            >> {log.out} 2>> {log.err}
        rm -f {input.r1} {input.r2}
        """

rule esviritu_mapped_to_assembly_database:
    input:
        r1="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_nonrrna.1.unmapped.fq.gz",
        r2="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}_nonrrna.2.unmapped.fq.gz",
        db_mmi="result/04_modified_esvirtue/databases/final_merged_database_only_assembly/virus_pathogen_database.mmi",
        db_metadata="result/04_modified_esvirtue/databases/final_merged_database_only_assembly/virus_pathogen_database.all_metadata.tsv",
    output:
        assembly_summary="result/04_modified_esvirtue/es/second_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv",
        info="result/04_modified_esvirtue/es/second_filter/{sample}/{sample}.detected_virus.info.tsv",
        consensus="result/04_modified_esvirtue/es/second_filter/{sample}/{sample}_final_consensus.fasta",
        coverage="result/04_modified_esvirtue/es/second_filter/{sample}/{sample}.virus_coverage_windows.tsv",
        log_file="result/04_modified_esvirtue/es/second_filter/{sample}/{sample}_esviritu.log",
        params="result/04_modified_esvirtue/es/second_filter/{sample}/{sample}_esviritu.params.yaml",
    conda:
        "../envs/esviritu_map.yaml"
    threads: config["esviritu_map"]["threads"]
    resources:
        mem_mb_per_cpu  = config["esviritu_map"]["memory"],
        runtime         = config["esviritu_map"]["runtime"],
        cpus_per_task   = config["esviritu_map"]["threads"],
        slurm_partition = config["esviritu_map"]["partition"],
        slurm_account   = config["esviritu_map"]["account"]
    params:
        quiet = config["esviritu_map"]["quiet"],
        extra = config["esviritu_map"].get("extra", "")
    log:
        out="log/04_modified_esvirtue/second_filter/{sample}.log",
        err="log/04_modified_esvirtue/second_filter/{sample}.err"
    shell:
        """
        mkdir -p result/04_modified_esvirtue/es/second_filter/{wildcards.sample}
        mkdir -p log/04_modified_esvirtue/second_filter
        EsViritu -r {input.r1} {input.r2} \
                 -s {wildcards.sample} \
                 -t {threads} \
                 -o result/04_modified_esvirtue/es/second_filter/{wildcards.sample} \
                 --db $(dirname {input.db_mmi}) -q {params.quiet} \
                 {params.extra} \
                 > {log.out} 2>> {log.err}
        rm -f {input.r1} {input.r2}
        rm -rf result/04_modified_esvirtue/es/first_filter/{wildcards.sample}/{wildcards.sample}_temp
        """

rule merge_esviritu_assembly_summary:
    input:
        first="result/04_modified_esvirtue/es/first_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv",
        second="result/04_modified_esvirtue/es/second_filter/{sample}/{sample}.detected_virus.assembly_summary.tsv",
    output:
        first_copy="result/04_modified_esvirtue/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.first.tsv",
        second_copy="result/04_modified_esvirtue/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.second.tsv",
        merged="result/04_modified_esvirtue/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.tsv",
    threads: config["merge_esviritu_assembly_summary"]["threads"]
    resources:
        mem_mb_per_cpu  = config["merge_esviritu_assembly_summary"]["memory"],
        runtime         = config["merge_esviritu_assembly_summary"]["runtime"],
        cpus_per_task   = config["merge_esviritu_assembly_summary"]["threads"],
        slurm_partition = config["merge_esviritu_assembly_summary"]["partition"],
        slurm_account   = config["merge_esviritu_assembly_summary"]["account"]
    log:
        out="log/04_modified_esvirtue/merge_info/{sample}.log",
        err="log/04_modified_esvirtue/merge_info/{sample}.err"
    shell:
        """
        mkdir -p result/04_modified_esvirtue/es/Merge/{wildcards.sample}
        mkdir -p log/04_modified_esvirtue/merge_info
        cp {input.first} {output.first_copy}
        cp {input.second} {output.second_copy}
        cat {input.first} {input.second} > {output.merged} 2> {log.err}
        """

rule rpkmf_esviritu_assembly_summary:
    input:
        merged="result/04_modified_esvirtue/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.tsv",
    output:
        rpkmf="result/04_modified_esvirtue/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.rpkmf.tsv",
    conda:
        "../envs/python.yaml"
    threads: config["rpkmf_esviritu_assembly_summary"]["threads"]
    resources:
        mem_mb_per_cpu  = config["rpkmf_esviritu_assembly_summary"]["memory"],
        runtime         = config["rpkmf_esviritu_assembly_summary"]["runtime"],
        cpus_per_task   = config["rpkmf_esviritu_assembly_summary"]["threads"],
        slurm_partition = config["rpkmf_esviritu_assembly_summary"]["partition"],
        slurm_account   = config["rpkmf_esviritu_assembly_summary"]["account"]
    params:
        script = config["scripts"]["esviritu_rpkmf"]
    log:
        out="log/04_modified_esvirtue/merge_info/{sample}.rpkmf.log",
        err="log/04_modified_esvirtue/merge_info/{sample}.rpkmf.err"
    shell:
        """
        python {params.script} --input {input.merged} --output {output.rpkmf} \
            > {log.out} 2> {log.err}
        """

rule merge_esviritu_rpkmf_all_samples:
    input:
        rpkmfs=expand(
            "result/04_modified_esvirtue/es/Merge/{sample}/{sample}.detected_virus.assembly_summary.merged.rpkmf.tsv",
            sample=config["rna_samples"],
        ),
    output:
        merged      = "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.tsv",
        subspecies  = "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.subspecies.tsv",
        read_count  = "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.read_count.tsv",
        subs_rc     = "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.subspecies.read_count.tsv",
        species_rc  = "result/04_modified_esvirtue/es/Merge/all_samples.detected_virus.assembly_summary.species.read_count.tsv",
    conda:
        "../envs/python.yaml"
    threads: config["merge_esviritu_rpkmf_all_samples"]["threads"]
    resources:
        mem_mb_per_cpu  = config["merge_esviritu_rpkmf_all_samples"]["memory"],
        runtime         = config["merge_esviritu_rpkmf_all_samples"]["runtime"],
        cpus_per_task   = config["merge_esviritu_rpkmf_all_samples"]["threads"],
        slurm_partition = config["merge_esviritu_rpkmf_all_samples"]["partition"],
        slurm_account   = config["merge_esviritu_rpkmf_all_samples"]["account"]
    params:
        script = config["scripts"]["esviritu_merge_rpkmf"]
    log:
        out="log/04_modified_esvirtue/merge_info/all_samples.rpkmf.log",
        err="log/04_modified_esvirtue/merge_info/all_samples.rpkmf.err"
    shell:
        """
        python {params.script} --inputs {input.rpkmfs} --output {output.merged} \
            > {log.out} 2> {log.err}
        """
