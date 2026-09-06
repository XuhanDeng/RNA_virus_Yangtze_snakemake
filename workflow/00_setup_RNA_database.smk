# RNA virus database setup workflow
# Run on login node (not via sbatch) — HPC compute nodes lack internet access.

configfile: "config/config.yaml"


rule all:
    input:
        directory(config["rdrp_catch"]["db_dir"]),
        directory(config["palm_annot"]["install_dir"]),
        directory(config["RVMT"]["dir"]),
        directory(config["esviritu_db"]["dir"]),
        # directory(config["kraken2_db"]["dir"]),
        # config["ICTV"]["fasta"],
        config["sortmerna"]["dbs"][0],
        config["sortmerna"]["dbs"][1],
        config["sortmerna"]["dbs"][2],
        config["sortmerna"]["dbs"][3],
        config["virushostdb"]["tsv"],
        config["checkv"]["db_marker"],
        config["genomad"]["db_marker"],
        config["virsorter2"]["db_marker"],
        config["ICTV_VMR"]["xlsx"],
        directory(config["diamond_nr"]["taxdump_dir"])
        # Note: generate_nr_fasta_and_taxid_map / build_diamond_nr_db moved to
        # workflow/00b_build_diamond_nr_db.smk (no internet needed there, so
        # it can run via sbatch on a compute node instead of the login node).





# Kraken2 database prepare



# ── SortMeRNA databases ───────────────────────────────────────────────────────

rule download_sortmerna_db:
    output:
        db1 = config["sortmerna"]["dbs"][0],
        db2 = config["sortmerna"]["dbs"][1],
        db3 = config["sortmerna"]["dbs"][2],
        db4 = config["sortmerna"]["dbs"][3]
    log:
        log = config["setup_log_dir"] + "/download_sortmerna_db/download_sortmerna_db.log",
        err = config["setup_log_dir"] + "/download_sortmerna_db/download_sortmerna_db.err"
    params:
        db_dir = config["sortmerna_db"]["dir"],
        url1   = config["sortmerna_db"]["smr_default_url"],
        url2   = config["sortmerna_db"]["smr_fast_url"],
        url3   = config["sortmerna_db"]["smr_sensitive_url"],
        url4   = config["sortmerna_db"]["smr_rfam_url"]
    shell:
        """
        mkdir -p {params.db_dir}
        wget -c '{params.url1}' -O {output.db1}.gz > {log.log} 2> {log.err}
        gunzip {output.db1}.gz >> {log.log} 2>> {log.err}
        wget -c '{params.url2}' -O {output.db2}.gz >> {log.log} 2>> {log.err}
        gunzip {output.db2}.gz >> {log.log} 2>> {log.err}
        wget -c '{params.url3}' -O {output.db3}.gz >> {log.log} 2>> {log.err}
        gunzip {output.db3}.gz >> {log.log} 2>> {log.err}
        wget -c '{params.url4}' -O {output.db4}.gz >> {log.log} 2>> {log.err}
        gunzip {output.db4}.gz >> {log.log} 2>> {log.err}
        """


# ── RdRpCATCH ─────────────────────────────────────────────────────────────────

rule download_rdrpcatch:
    output:
        directory(config["rdrp_catch"]["db_dir"])
    log:
        log = config["setup_log_dir"] + "/download_rdrpcatch/download_rdrpcatch.log",
        err = config["setup_log_dir"] + "/download_rdrpcatch/download_rdrpcatch.err"
    params:
        db_dir = config["rdrp_catch"]["db_dir"]
    conda:
        "../envs/rdrp_catch.yaml"
    shell:
        """
        rdrpcatch databases --destination-dir {params.db_dir} \
            > {log.log} 2> {log.err}
        """


# ── palm_annot ────────────────────────────────────────────────────────────────
# Requires Ubuntu/Linux — precompiled binaries, no conda env needed.
# Github: https://github.com/rcedgar/palm_annot

rule install_palm_annot:
    output:
        directory(config["palm_annot_install"]["install_dir"])
    log:
        log = config["setup_log_dir"] + "/install_palm_annot/install_palm_annot.log",
        err = config["setup_log_dir"] + "/install_palm_annot/install_palm_annot.err"
    params:
        url         = config["palm_annot_install"]["url"],
        install_dir = config["palm_annot_install"]["install_dir"]
    shell:
        """
        git clone {params.url} {params.install_dir} \
            > {log.log} 2> {log.err}
        chmod +x {params.install_dir}/bin/* {params.install_dir}/py/* \
            >> {log.log} 2>> {log.err}
        """


# ── RVMT ──────────────────────────────────────────────────────────────────────

rule download_RVMT:
    output:
        config["RVMT"]["zip"]
    log:
        log = config["setup_log_dir"] + "/download_RVMT/download_RVMT.log",
        err = config["setup_log_dir"] + "/download_RVMT/download_RVMT.err"
    params:
        url = config["RVMT"]["url"]
    shell:
        "wget -c '{params.url}' -O {output} > {log.log} 2> {log.err}"


rule extract_RVMT:
    input:
        config["RVMT"]["zip"]
    output:
        directory(config["RVMT"]["dir"])
    log:
        log = config["setup_log_dir"] + "/extract_RVMT/extract_RVMT.log",
        err = config["setup_log_dir"] + "/extract_RVMT/extract_RVMT.err"
    shell:
        "unzip {input} -d {output} > {log.log} 2> {log.err}"


# ── ICTV ──────────────────────────────────────────────────────────────────────

# rule download_ICTV_vmr:
#     output:
#         config["ICTV"]["vmr_xlsx"]
#     log:
#         log = config["setup_log_dir"] + "/download_ICTV_vmr/download_ICTV_vmr.log",
#         err = config["setup_log_dir"] + "/download_ICTV_vmr/download_ICTV_vmr.err"
#     params:
#         url = config["ICTV"]["vmr_url"]
#     shell:
#         "wget -L -c '{params.url}' -O {output} > {log.log} 2> {log.err}"
#
#
# rule extract_ICTV_accessions:
#     input:
#         config["ICTV"]["vmr_xlsx"]
#     output:
#         config["ICTV"]["accessions_txt"]
#     log:
#         log = config["setup_log_dir"] + "/extract_ICTV_accessions/extract_ICTV_accessions.log",
#         err = config["setup_log_dir"] + "/extract_ICTV_accessions/extract_ICTV_accessions.err"
#     params:
#         sheet           = config["ICTV"]["vmr_sheet"],
#         realm_filter    = config["ICTV"]["realm_filter"],
#         coverage_filter = config["ICTV"]["coverage_filter"]
#     conda:
#         "../envs/openpyxl.yaml"
#     script:
#         "../scripts/00_setup_RNA_database/extract_ICTV_accessions.py"
#
#
# rule download_ICTV_sequences:
#     input:
#         config["ICTV"]["accessions_txt"]
#     output:
#         config["ICTV"]["fasta"]
#     log:
#         log = config["setup_log_dir"] + "/download_ICTV_sequences/download_ICTV_sequences.log",
#         err = config["setup_log_dir"] + "/download_ICTV_sequences/download_ICTV_sequences.err"
#     conda:
#         "../envs/entrez_direct.yaml"
#     script:
#         "../scripts/00_setup_RNA_database/download_ICTV_sequences.py"


# ── ICTV VMR (for GeNomad virus annotation) ──────────────────────────────────

rule download_ICTV_VMR:
    output:
        config["ICTV_VMR"]["xlsx"]
    log:
        err = config["setup_log_dir"] + "/download_ICTV_VMR.err"
    params:
        url = config["ICTV_VMR"]["url"]
    shell:
        """
        mkdir -p $(dirname {output})
        wget -c '{params.url}' -O {output} 2> {log.err}
        """


# ── ESViritu database ─────────────────────────────────────────────────────────

rule download_esviritu_db:
    output:
        directory(config["esviritu_db"]["dir"])
    log:
        log = config["setup_log_dir"] + "/download_esviritu_db/download_esviritu_db.log",
        err = config["setup_log_dir"] + "/download_esviritu_db/download_esviritu_db.err"
    params:
        url    = config["esviritu_db"]["url"],
        tar    = config["esviritu_db"]["tar"],
        db_dir = config["esviritu_db"]["dir"]
    shell:
        """
        mkdir -p {params.db_dir}
        wget -c '{params.url}' -O {params.tar} > {log.log} 2> {log.err}
        tar -xzf {params.tar} -C {params.db_dir} --strip-components=1 >> {log.log} 2>> {log.err}
        """


# ── Kraken2 database ──────────────────────────────────────────────────────────

rule download_kraken2_db:
    output:
        directory(config["kraken2_db"]["dir"])
    log:
        log = config["setup_log_dir"] + "/download_kraken2_db/download_kraken2_db.log",
        err = config["setup_log_dir"] + "/download_kraken2_db/download_kraken2_db.err"
    params:
        url    = config["kraken2_db"]["url"],
        tar    = config["kraken2_db"]["tar"],
        db_dir = config["kraken2_db"]["dir"]
    shell:
        """
        mkdir -p {params.db_dir}
        wget -c '{params.url}' -O {params.tar} > {log.log} 2> {log.err}
        tar -xzf {params.tar} -C {params.db_dir} >> {log.log} 2>> {log.err}
        """


# ── Virus-Pathogen DB host metadata ──────────────────────────────────────────

rule download_virus_pathogen_db_host_metadata:
    output:
        zip = config["virus_pathogen_db"]["zip"],
        tsv = config["virus_pathogen_db"]["tsv"]
    log:
        log = config["setup_log_dir"] + "/download_virus_pathogen_db_host_metadata/download.log",
        err = config["setup_log_dir"] + "/download_virus_pathogen_db_host_metadata/download.err"
    params:
        url = config["virus_pathogen_db"]["url"],
        dir = lambda wildcards, output: os.path.dirname(output.zip)
    shell:
        """
        mkdir -p {params.dir}
        wget -c '{params.url}' -O {output.zip} > {log.log} 2> {log.err}
        unzip -o {output.zip} -d {params.dir} >> {log.log} 2>> {log.err}
        """


# ── Virus-Host DB ─────────────────────────────────────────────────────────────

rule download_virushostdb:
    output:
        config["virushostdb"]["tsv"]
    log:
        log = config["setup_log_dir"] + "/download_virushostdb/download_virushostdb.log",
        err = config["setup_log_dir"] + "/download_virushostdb/download_virushostdb.err"
    params:
        url = config["virushostdb"]["url"],
        dir = lambda wildcards, output: os.path.dirname(output[0])
    shell:
        """
        mkdir -p {params.dir}
        wget -c '{params.url}' -O {output} > {log.log} 2> {log.err}
        """





## ssDNA realated


# ── ssDNA-related databases ───────────────────────────────────────────────────

# Database 1: CheckV
rule download_checkv_db:
    output:
        marker = config["checkv"]["db_marker"]
    conda:
        "../envs/checkv.yaml"
    log:
        log = config["setup_log_dir"] + "/download_checkv_db/download_checkv_db.log",
        err = config["setup_log_dir"] + "/download_checkv_db/download_checkv_db.err"
    threads: config["checkv"]["threads"]
    params:
        db_dir = config["checkv"]["db_dir"]
    shell:
        """
        mkdir -p {params.db_dir}
        checkv download_database {params.db_dir} > {log.log} 2> {log.err}
        touch {output.marker}
        """


# Database 2: geNomad
rule download_genomad_db:
    output:
        marker = config["genomad"]["db_marker"]
    conda:
        "../envs/genomad.yaml"
    log:
        log = config["setup_log_dir"] + "/download_genomad_db/download_genomad_db.log",
        err = config["setup_log_dir"] + "/download_genomad_db/download_genomad_db.err"
    threads: config["genomad"]["threads"]
    params:
        download_dir = config["genomad"]["download_dir"]
    shell:
        """
        mkdir -p {params.download_dir}
        genomad download-database {params.download_dir} > {log.log} 2> {log.err}
        touch {output.marker}
        """


# Database 3: VirSorter2
rule download_virsorter2_db:
    output:
        marker = config["virsorter2"]["db_marker"]
    conda:
        "../envs/vs2.yaml"
    log:
        log = config["setup_log_dir"] + "/download_virsorter2_db/download_virsorter2_db.log",
        err = config["setup_log_dir"] + "/download_virsorter2_db/download_virsorter2_db.err"
    threads: config["virsorter2"]["threads"]
    params:
        db_dir = config["virsorter2"]["database"]
    shell:
        """
        mkdir -p {params.db_dir}
        virsorter setup -d {params.db_dir} -j {threads} > {log.log} 2> {log.err}
        touch {output.marker}
        """


# ── Diamond nr database (viral-origin check, workflow 03 Step 10) ───────────
# NCBI stopped publishing the plain nr.gz FASTA on the FTP site in April
# 2024 (https://ncbiinsights.ncbi.nlm.nih.gov/2024/01/25/blast-fasta-unavailable-on-ftp/).
# The only supported source now is the pre-formatted BLAST database, fetched
# via update_blastdb.pl (ships with blast+) and regenerated to FASTA (+ an
# accession->taxid mapping) with blastdbcmd -- both come from the same BLAST
# db volumes, so no separate prot.accession2taxid.FULL.gz download needed.
# The resulting FASTA + mapping + NCBI taxdump are used to build a Diamond
# database with --taxonmap/--taxonnodes/--taxonnames, which lets diamond
# blastp filter searches directly by taxon at query time (--taxonlist
# <taxid>) and report per-hit staxids -- no separate taxonomy library (e.g.
# ete3) needed anywhere downstream.
#
# NOTE: nr is very large (100+ GB, growing). Run on a node/partition with
# enough disk and RAM, and expect this to take hours.

rule download_blast_nr_db:
    output:
        marker = config["diamond_nr"]["db_dir"] + "/.blastdb_download.done"
    log:
        log = config["setup_log_dir"] + "/download_blast_nr_db/download.log",
        err = config["setup_log_dir"] + "/download_blast_nr_db/download.err"
    conda:
        "../envs/blast.yaml"
    threads: config["diamond_nr"]["threads"]
    resources:
        mem_mb_per_cpu  = config["diamond_nr"]["memory"],
        runtime         = config["diamond_nr"]["runtime"],
        cpus_per_task   = config["diamond_nr"]["threads"],
        slurm_partition = config["diamond_nr"]["partition"],
        slurm_account   = config["diamond_nr"]["account"]
    params:
        db_dir = config["diamond_nr"]["db_dir"]
    shell:
        """
        mkdir -p {params.db_dir}
        mkdir -p $(dirname {log.log})
        log_abs=$(realpath {log.log})
        err_abs=$(realpath {log.err})
        marker_abs=$(realpath -m {output.marker})
        cd {params.db_dir}
        update_blastdb.pl --decompress --num_threads {threads} nr \
            > $log_abs 2> $err_abs
        touch $marker_abs
        """


rule download_diamond_nr_taxdump:
    output:
        directory(config["diamond_nr"]["taxdump_dir"])
    log:
        log = config["setup_log_dir"] + "/download_diamond_nr_taxdump/download.log",
        err = config["setup_log_dir"] + "/download_diamond_nr_taxdump/download.err"
    params:
        url = config["diamond_nr"]["taxdump_url"],
        dir = config["diamond_nr"]["taxdump_dir"]
    shell:
        """
        mkdir -p {params.dir}
        mkdir -p $(dirname {log.log})
        wget -c '{params.url}' -O {params.dir}/new_taxdump.tar.gz > {log.log} 2> {log.err}
        tar -xzf {params.dir}/new_taxdump.tar.gz -C {params.dir} >> {log.log} 2>> {log.err}
        """
