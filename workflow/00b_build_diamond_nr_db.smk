# Diamond nr database build (no internet required)
# Run via sbatch on a compute node, AFTER workflow/00_setup_RNA_database.smk
# has finished download_blast_nr_db + download_diamond_nr_taxdump on the
# login node (both of those need internet; these two rules only read
# already-downloaded local files, so they belong on a compute node instead).

configfile: "config/config.yaml"


rule all:
    input:
        config["diamond_nr"]["db_marker"],
        config["diamond_nr"]["db_dir"] + "/.blastdb_cleanup.done"


# Step 1a: extract FASTA from BLAST db (single-threaded, I/O bound, ~12-24h)
rule generate_nr_fasta:
    input:
        marker = config["diamond_nr"]["db_dir"] + "/.blastdb_download.done"
    output:
        nr_fasta = config["diamond_nr"]["nr_fasta"]
    log:
        log = config["setup_log_dir"] + "/generate_nr_fasta/generate.log",
        err = config["setup_log_dir"] + "/generate_nr_fasta/generate.err"
    conda:
        "../envs/blast.yaml"
    threads: 2
    resources:
        mem_mb_per_cpu  = config["diamond_nr"]["blastdbcmd_memory"],
        runtime         = config["diamond_nr"]["runtime"],
        cpus_per_task   = 2,
        slurm_partition = config["diamond_nr"]["blastdbcmd_partition"],
        slurm_account   = config["diamond_nr"]["account"]
    params:
        db_dir = config["diamond_nr"]["db_dir"]
    shell:
        """
        mkdir -p $(dirname {log.log}) $(dirname {output.nr_fasta})
        cd {params.db_dir}
        blastdbcmd -db nr -entry all -out $(realpath -m {output.nr_fasta}) \
            > $(realpath {log.log}) 2> $(realpath {log.err})
        """


# Step 1b: extract accession->taxid map (single-threaded, I/O bound, ~12-24h)
# Runs in parallel with generate_nr_fasta as a separate SLURM job
rule generate_nr_taxid_map:
    input:
        marker = config["diamond_nr"]["db_dir"] + "/.blastdb_download.done"
    output:
        acc2tax = config["diamond_nr"]["prot_accession2taxid"]
    log:
        err = config["setup_log_dir"] + "/generate_nr_taxid_map/generate.err"
    conda:
        "../envs/blast.yaml"
    threads: 2
    resources:
        mem_mb_per_cpu  = config["diamond_nr"]["blastdbcmd_memory"],
        runtime         = config["diamond_nr"]["runtime"],
        cpus_per_task   = 2,
        slurm_partition = config["diamond_nr"]["blastdbcmd_partition"],
        slurm_account   = config["diamond_nr"]["account"]
    params:
        db_dir = config["diamond_nr"]["db_dir"]
    shell:
        """
        mkdir -p $(dirname {log.err}) $(dirname {output.acc2tax})
        cd {params.db_dir}
        echo -e "accession.version\ttaxid" > $(realpath -m {output.acc2tax})
        blastdbcmd -db nr -entry all -outfmt "%a %T" \
            | awk '{{print $1"\\t"$2}}' >> $(realpath -m {output.acc2tax}) \
            2> $(realpath {log.err})
        """


# Step 2: build diamond database (multi-threaded, ~2-4h)
rule build_diamond_nr_db:
    input:
        nr_fasta = config["diamond_nr"]["nr_fasta"],
    output:
        marker = config["diamond_nr"]["db_marker"]
    conda:
        "../envs/diamond.yaml"
    log:
        log = config["setup_log_dir"] + "/build_diamond_nr_db/build.log",
        err = config["setup_log_dir"] + "/build_diamond_nr_db/build.err"
    threads: config["diamond_nr"]["threads"]
    resources:
        mem_mb_per_cpu  = config["diamond_nr"]["memory"],
        runtime         = config["diamond_nr"]["runtime"],
        cpus_per_task   = config["diamond_nr"]["threads"],
        slurm_partition = config["diamond_nr"]["partition"],
        slurm_account   = config["diamond_nr"]["account"]
    params:
        dmnd_db = config["diamond_nr"]["dmnd_db"]
    shell:
        """
        mkdir -p $(dirname {params.dmnd_db})
        diamond makedb \
            --in {input.nr_fasta} \
            -d {params.dmnd_db} \
            -p {threads} \
            > {log.log} 2> {log.err}
        touch {output.marker}
        """


# Step 3: delete BLAST db volumes, tarballs, and nr.fasta (~1.7 TB freed)
rule cleanup_nr_blastdb:
    input:
        dmnd_marker = config["diamond_nr"]["db_marker"],
    output:
        marker = config["diamond_nr"]["db_dir"] + "/.blastdb_cleanup.done"
    log:
        err = config["setup_log_dir"] + "/cleanup_nr_blastdb/cleanup.err"
    threads: 1
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = 1,
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"]
    params:
        db_dir = config["diamond_nr"]["db_dir"]
    shell:
        """
        mkdir -p $(dirname {log.err})
        touch {log.err}
        ABS_LOG=$(realpath {log.err})
        ABS_MARKER=$(realpath -m {output.marker})
        cd {params.db_dir}
        rm -f nr.fasta \
              nr.*.tar.gz nr.*.tar.gz.md5 \
              nr.*.phr nr.*.psq nr.*.pin nr.*.pog nr.*.psi nr.*.psd \
              nr.pjs nr.pot nr.ptf nr.pto nr.pal nr.pdb nr.pos \
              taxdb.btd taxdb.bti taxonomy4blast.sqlite3 \
              2> "$ABS_LOG"
        touch "$ABS_MARKER"
        """
