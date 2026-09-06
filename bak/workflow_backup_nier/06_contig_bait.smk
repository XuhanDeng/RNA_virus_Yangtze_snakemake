configfile: "config/config.yaml"

# ------------------------------------------------------------------ #
# Bait-based viral contig recovery                                    #
#                                                                     #
# Query (bait): RVMT NovoContigs + tier1a/1b/2 RdRP confirmed contigs #
# Target:       all assembled contigs from workflow 02 (all samples   #
#               cat-ed into one FASTA)                                #
#                                                                     #
# Method: mmseqs search (nucleotide, non-sensitive) followed by       #
# stringent identity/coverage filtering to recover highly similar     #
# contigs from the bulk set.                                          #
# ------------------------------------------------------------------ #

_OUTDIR      = "result/06_contig_bait"
_BAIT_DIR    = f"{_OUTDIR}/1_bait_db"
_TARGET_DIR  = f"{_OUTDIR}/2_target_db"
_SEARCH_DIR  = f"{_OUTDIR}/3_mmseqs_search"
_RESULT_DIR  = f"{_OUTDIR}/4_filtered_hits"
_CONTIG_DIR  = f"{_OUTDIR}/5_bait_contigs"
_TMP_DIR     = f"{_OUTDIR}/tmp"

_NOVO_CONTIGS   = "database/RVMT/RVMT_Zenodo_V4/RVMT_Zenodo_V4/Expanded_Contig_Set/NovoContigs.fasta"
_MOTIF_RESULTS   = "result/03_RDRP_identification/4_motif_search/motif_results"
_TIER_SUMMARY    = _MOTIF_RESULTS + "/tier_summary.tsv"
_BAIT_RDRP_FASTA  = f"{_OUTDIR}/0_all_contigs/bait_rdrp_contigs.fasta"
_TARGET_CONTIGS   = f"{_OUTDIR}/0_all_contigs/target_contigs.fasta"

# All assembled contigs from workflow 02 (one FASTA per sample, cat-ed here)
_RENAME_DIR     = config["rna_reformated_scaffolds_dir"]
_MIN_LEN        = config["seqkit"]["min_length"]
_ALL_CONTIGS    = f"{_OUTDIR}/0_all_contigs/all_samples.fasta"
_ALL_CONTIG_FASTAS = expand(
    _RENAME_DIR + "/rename_{min_len}/{sample}_scaffolds_rename_{min_len}.fasta",
    sample=config["rna_samples"],
    min_len=_MIN_LEN,
)


rule all:
    input:
        _ALL_CONTIGS,
        f"{_CONTIG_DIR}/bait_recovered.txt",
        f"{_CONTIG_DIR}/bait_recovered.fasta",
        f"{_RESULT_DIR}/hits.filtered.tsv",
        f"{_RESULT_DIR}/bait_target_map.tsv",
        # hits.filtered.tsv  — all passing alignments (one row per hit)
        # bait_target_map.tsv — one row per recovered contig: best bait hit per source


# ------------------------------------------------------------------ #
# Step 0: Combine all workflow 02 contigs into one target FASTA       #
# ------------------------------------------------------------------ #

rule combine_all_contigs:
    input:
        ancient(_ALL_CONTIG_FASTAS),
    output:
        fasta = _ALL_CONTIGS,
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        err = "log/06_contig_bait/0_combine_contigs.err",
    shell:
        """
        mkdir -p $(dirname {output.fasta}) log/06_contig_bait
        cat {input} > {output.fasta} 2> {log.err}
        echo "Total contigs: $(grep -c '^>' {output.fasta})" >> {log.err}
        """


# ------------------------------------------------------------------ #
# Step 0b: Extract tier1a/1b/2 nucleotide contigs from all-contigs   #
#          using contig IDs from tier_summary.tsv                     #
# ------------------------------------------------------------------ #

rule extract_rdrp_bait_contigs:
    input:
        tsv   = ancient(_TIER_SUMMARY),
        fasta = _ALL_CONTIGS,
    output:
        bait   = _BAIT_RDRP_FASTA,
        target = _TARGET_CONTIGS,
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        out = "log/06_contig_bait/0b_extract_rdrp_bait.log",
        err = "log/06_contig_bait/0b_extract_rdrp_bait.err",
    shell:
        """
        # extract tier1a/1b/2 contig IDs (col2=tier, col3=contig)
        awk -F'\\t' 'NR>1 && ($2=="tier1a" || $2=="tier1b" || $2=="tier2") {{print $3}}' \
            {input.tsv} | sort -u > {output.bait}.ids.tmp
        # bait: our RdRp contigs (nucleotide)
        seqkit grep \
            --pattern-file {output.bait}.ids.tmp \
            --threads {threads} \
            {input.fasta} \
            > {output.bait} \
            2> {log.err}
        # target: all contigs excluding bait IDs
        seqkit grep \
            --pattern-file {output.bait}.ids.tmp \
            --invert-match \
            --threads {threads} \
            {input.fasta} \
            > {output.target} \
            2>> {log.err}
        echo "RdRp bait contigs: $(grep -c '^>' {output.bait})" > {log.out}
        echo "Target contigs:    $(grep -c '^>' {output.target})" >> {log.out}
        rm -f {output.bait}.ids.tmp
        """


# ------------------------------------------------------------------ #
# Step 1: Combine bait sequences (RVMT NovoContigs + tier1a/1b/2)    #
# ------------------------------------------------------------------ #

rule combine_bait:
    input:
        novo = ancient(_NOVO_CONTIGS),
        rdrp = _BAIT_RDRP_FASTA,
    output:
        combined = f"{_BAIT_DIR}/bait_combined.fasta",
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        out = "log/06_contig_bait/1_combine_bait.log",
        err = "log/06_contig_bait/1_combine_bait.err",
    shell:
        """
        mkdir -p {_BAIT_DIR} log/06_contig_bait
        cat {input.novo} {input.rdrp} > {output.combined} 2> {log.err}
        echo "Combined bait fasta written to {output.combined}" > {log.out}
        """


# ------------------------------------------------------------------ #
# Step 2: Build mmseqs2 databases                                    #
# ------------------------------------------------------------------ #

rule mmseqs_createdb_query:
    input:
        fasta = f"{_BAIT_DIR}/bait_combined.fasta",
    output:
        db = f"{_BAIT_DIR}/bait_db.index",
    conda:
        "../envs/mmseqs2.yaml"
    threads: config["contig_bait"]["threads"]
    resources:
        mem_mb_per_cpu  = config["contig_bait"]["memory"],
        runtime         = config["contig_bait"]["runtime"],
        cpus_per_task   = config["contig_bait"]["threads"],
        slurm_partition = config["contig_bait"]["partition"],
        slurm_account   = config["contig_bait"]["account"],
    log:
        out = "log/06_contig_bait/2_createdb_query.log",
        err = "log/06_contig_bait/2_createdb_query.err",
    shell:
        """
        mmseqs createdb {input.fasta} {_BAIT_DIR}/bait_db > {log.out} 2> {log.err}
        """


rule mmseqs_createdb_target:
    input:
        fasta = _TARGET_CONTIGS,
    output:
        db = f"{_TARGET_DIR}/target_db.index",
    conda:
        "../envs/mmseqs2.yaml"
    threads: config["contig_bait"]["threads"]
    resources:
        mem_mb_per_cpu  = config["contig_bait"]["memory"],
        runtime         = config["contig_bait"]["runtime"],
        cpus_per_task   = config["contig_bait"]["threads"],
        slurm_partition = config["contig_bait"]["partition"],
        slurm_account   = config["contig_bait"]["account"],
    log:
        out = "log/06_contig_bait/2_createdb_target.log",
        err = "log/06_contig_bait/2_createdb_target.err",
    shell:
        """
        mkdir -p {_TARGET_DIR}
        mmseqs createdb {input.fasta} {_TARGET_DIR}/target_db > {log.out} 2> {log.err}
        """


# ------------------------------------------------------------------ #
# Step 3: mmseqs search (nucleotide, non-sensitive)                  #
# Parameters from paper:                                             #
#   --search-type 3  nucleotide-nucleotide                           #
#   --min-aln-len 120                                                #
#   --min-seq-id 0.66  (loose initial filter; stringent post-filter) #
#   -s 1               non-sensitive                                 #
#   -c 0.85 --cov-mode 1  target coverage                           #
# ------------------------------------------------------------------ #

rule mmseqs_search:
    input:
        query_db  = f"{_BAIT_DIR}/bait_db.index",
        target_db = f"{_TARGET_DIR}/target_db.index",
    output:
        result_db = f"{_SEARCH_DIR}/result_db.index",
    conda:
        "../envs/mmseqs2.yaml"
    threads: config["contig_bait"]["threads"]
    resources:
        mem_mb_per_cpu  = config["contig_bait"]["memory"],
        runtime         = config["contig_bait"]["runtime"],
        cpus_per_task   = config["contig_bait"]["threads"],
        slurm_partition = config["contig_bait"]["partition"],
        slurm_account   = config["contig_bait"]["account"],
    log:
        out = "log/06_contig_bait/3_mmseqs_search.log",
        err = "log/06_contig_bait/3_mmseqs_search.err",
    shell:
        """
        mkdir -p {_SEARCH_DIR} {_TMP_DIR}
        mmseqs search \
            {_BAIT_DIR}/bait_db \
            {_TARGET_DIR}/target_db \
            {_SEARCH_DIR}/result_db \
            {_TMP_DIR} \
            --search-type 3 \
            --min-aln-len 120 \
            --min-seq-id 0.66 \
            -s 1 \
            -c 0.85 \
            --cov-mode 1 \
            --threads {threads} \
            > {log.out} 2> {log.err}
        """


# ------------------------------------------------------------------ #
# Step 4: Convert result to readable TSV                             #
# ------------------------------------------------------------------ #

rule mmseqs_convertalis:
    input:
        query_db  = f"{_BAIT_DIR}/bait_db.index",
        target_db = f"{_TARGET_DIR}/target_db.index",
        result_db = f"{_SEARCH_DIR}/result_db.index",
    output:
        tsv = f"{_SEARCH_DIR}/hits.tsv",
    conda:
        "../envs/mmseqs2.yaml"
    threads: config["contig_bait"]["threads"]
    resources:
        mem_mb_per_cpu  = config["contig_bait"]["memory"],
        runtime         = config["contig_bait"]["runtime"],
        cpus_per_task   = config["contig_bait"]["threads"],
        slurm_partition = config["contig_bait"]["partition"],
        slurm_account   = config["contig_bait"]["account"],
    log:
        out = "log/06_contig_bait/4_convertalis.log",
        err = "log/06_contig_bait/4_convertalis.err",
    shell:
        """
        mmseqs convertalis \
            {_BAIT_DIR}/bait_db \
            {_TARGET_DIR}/target_db \
            {_SEARCH_DIR}/result_db \
            {output.tsv} \
            --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qlen,tlen,tcov" \
            --threads {threads} \
            > {log.out} 2> {log.err}
        """


# ------------------------------------------------------------------ #
# Step 5: Stringent post-filter                                      #
#   E-value < 1e-9, Identity > 95%, target coverage >= 95%          #
# ------------------------------------------------------------------ #

rule filter_bait_hits:
    input:
        tsv          = f"{_SEARCH_DIR}/hits.tsv",
        tier_summary = ancient(_TIER_SUMMARY),
    output:
        filtered = f"{_RESULT_DIR}/hits.filtered.tsv",
        bait_map = f"{_RESULT_DIR}/bait_target_map.tsv",
    conda:
        "../envs/python.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    params:
        script      = "scripts/06_contig_bait/filter_bait_hits.py",
        evalue      = config["contig_bait"]["filter_evalue"],
        min_pident  = config["contig_bait"]["filter_min_pident"],
        min_tcov    = config["contig_bait"]["filter_min_tcov"],
    log:
        out = "log/06_contig_bait/5_filter_hits.log",
        err = "log/06_contig_bait/5_filter_hits.err",
    shell:
        """
        mkdir -p {_RESULT_DIR}
        python {params.script} \
            --input          {input.tsv} \
            --output         {output.filtered} \
            --bait-map       {output.bait_map} \
            --tier-summary   {input.tier_summary} \
            --evalue         {params.evalue} \
            --min-pident     {params.min_pident} \
            --min-tcov       {params.min_tcov} \
            > {log.out} 2> {log.err}
        """


# ------------------------------------------------------------------ #
# Step 6: Extract recovered contigs from no_rdrp.fasta               #
# ------------------------------------------------------------------ #

rule extract_bait_contigs:
    input:
        filtered = f"{_RESULT_DIR}/hits.filtered.tsv",
        fasta    = _TARGET_CONTIGS,
    output:
        txt   = f"{_CONTIG_DIR}/bait_recovered.txt",
        fasta = f"{_CONTIG_DIR}/bait_recovered.fasta",
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        mem_mb_per_cpu  = config["seqkit"]["memory"],
        runtime         = config["seqkit"]["runtime"],
        cpus_per_task   = config["seqkit"]["threads"],
        slurm_partition = config["seqkit"]["partition"],
        slurm_account   = config["seqkit"]["account"],
    log:
        out = "log/06_contig_bait/6_extract_contigs.log",
        err = "log/06_contig_bait/6_extract_contigs.err",
    shell:
        """
        mkdir -p {_CONTIG_DIR}
        awk 'NR>1 {{print $2}}' {input.filtered} | sort -u > {output.txt} 2> {log.err}
        seqkit grep \
            --pattern-file {output.txt} \
            --threads {threads} \
            {input.fasta} \
            > {output.fasta} \
            2>> {log.err}
        echo "Recovered $(wc -l < {output.txt}) contigs" > {log.out}
        """
