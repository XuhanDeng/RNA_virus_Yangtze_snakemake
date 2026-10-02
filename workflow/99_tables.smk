# Summary table generation workflow
# Collects key outputs from upstream workflows and produces final publication tables.
# Run after all upstream workflows are complete.
#
# Trimmed down from the original 99_tables.smk (full version backed up at
# 99_tables.smk.bak) -- RNA virus contig table, recalculated RPKM/correlation,
# read stats, and SparCC all removed for now. Bait recovery is being
# abandoned, and the RdRp-based recalculation will be rebuilt separately with
# its own spec later.

configfile: "config/config.yaml"

# ── Input → output mapping ────────────────────────────────────────────────────

_TABLES = {
    # ESvirtu abundance
    "04_modified_esvirtue_ribodetector/esvirtu_read_count.tsv":
        "result/04_modified_esvirtue_ribodetector/es/Merge/all_samples.detected_virus.assembly_summary.read_count.tsv",
    "04_modified_esvirtue_ribodetector/esvirtu_read_count.species.tsv":
        "result/04_modified_esvirtue_ribodetector/es/Merge/all_samples.detected_virus.assembly_summary.species.read_count.tsv",
    "04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv":
        "result/04_modified_esvirtue_ribodetector/es/Merge/all_samples.detected_virus.assembly_summary.subspecies.read_count.tsv",
    "04_modified_esvirtue_ribodetector/esvirtu_rpkmf.tsv":
        "result/04_modified_esvirtue_ribodetector/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.tsv",
    "04_modified_esvirtue_ribodetector/esvirtu_rpkmf.species.tsv":
        "result/04_modified_esvirtue_ribodetector/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.species.tsv",
    "04_modified_esvirtue_ribodetector/esvirtu_rpkmf.subspecies.tsv":
        "result/04_modified_esvirtue_ribodetector/es/Merge/all_samples.detected_virus.assembly_summary.rpkmf.subspecies.tsv",
    # ESvirtu reference
    "03b_esvirtu_rdrp_identification/esvirtu_reference_metadata.tsv":
        "result/03b_esvirtu_rdrp_identification/2_metadata/esvirtu_reference.metadata.tsv",
    # Other virus identification (GeNomad + ICTV VMR annotation)
    "08_other_virus_identification/all_samples_virus_filter_summary.annotated.tsv":
        "result/08_other_virus_identification/2_virus_filter_summary/all_samples_virus_filter_summary.annotated.tsv",
    # RVMT DIAMOND taxonomy on the CD-HIT cluster representatives (Step 21) --
    # destination mirrors the source path exactly (21_taxonomic_assignment_cluster/full_length/)
    "03_RDRP_identification/21_taxonomic_assignment_cluster/full_length/full_length_diamond_rvmt_annotated.tsv":
        "result/03_RDRP_identification/21_taxonomic_assignment_cluster/full_length/full_length_diamond_rvmt_annotated.tsv",
    "03_RDRP_identification/21_taxonomic_assignment_cluster/full_length/full_length_diamond_rvmt_tophit.tsv":
        "result/03_RDRP_identification/21_taxonomic_assignment_cluster/full_length/full_length_diamond_rvmt_tophit.tsv",
    # Final protein summary (Step 8)
    "03_RDRP_identification/8_RdRp_protein/final_proteins_summary.tsv":
        "result/03_RDRP_identification/8_RdRp_protein/final_proteins_summary.tsv",
    # CD-HIT cluster membership file (Step 20) -- note: the actual result dir
    # is 10_final/nr_filtered_cluster/, not 20_cluster/ (20_cluster only
    # appears as a log path in 03_RDRP_identification.smk).
    "03_RDRP_identification/10_final/nr_filtered_cluster/combined_full_c90.faa.clstr":
        "result/03_RDRP_identification/10_final/nr_filtered_cluster/combined_full_c90.faa.clstr",
    # Final NR-confirmed per-contig RdRp identification table (Step 10) -- one
    # row per contig, source (RdRpCATCH/LucaProt), motif region coords, etc.
    "03_RDRP_identification/10_final/final_rdrp_merged.tsv":
        "result/03_RDRP_identification/10_final/final_rdrp_merged.tsv",
}

_OUTDIR  = "result/99_tables/01_cp_table"
_LOG_DIR = "log/99_tables"

_ESVIRTU_INFO_TABLE = "result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv"
_RNA_VIRUS_TABLE    = "result/99_tables/4_contig_tag_abundance/rna_virus_contig_table.tsv"

# ── Recalculated RNA-virus-only abundance (esvirtu ssRNA/dsRNA rows +
# RdRpCATCH/LucaProt contig rows) ──────────────────────────────────────────────
_RECALC_OUTDIR = "result/99_tables/5_Recalculated_RPKM_result/0_RDRP_Esvirtue"
_RECALC_COUNT  = _RECALC_OUTDIR + "/rna_virus_recalc_read_count.tsv"
_RECALC_RPKMF  = _RECALC_OUTDIR + "/rna_virus_recalc_rpkmf.tsv"
_RECALC_TPM    = _RECALC_OUTDIR + "/rna_virus_recalc_tpm.tsv"

# ── RVMT taxonomy annotation of the recalculated tables (contig rows only;
# esvirtu rows keep their own pre-existing taxonomy, untouched) ───────────────
_TAXONOMY_OUTDIR = "result/99_tables/5_Recalculated_RPKM_result/1_taxonmy"
_TAX_COUNT = _TAXONOMY_OUTDIR + "/rna_virus_recalc_read_count.annotated.tsv"
_TAX_RPKMF = _TAXONOMY_OUTDIR + "/rna_virus_recalc_rpkmf.annotated.tsv"
_TAX_TPM   = _TAXONOMY_OUTDIR + "/rna_virus_recalc_tpm.annotated.tsv"

# ── Per-sequence length/GC stats for the 03_RDRP_identification 10_final
# all/nr_filtered/nr_filtered_cluster contig + protein outputs ────────────────
_RDRP_STAGE_FASTA = {
    ("all",                "rdrp_contigs"): "result/03_RDRP_identification/10_final/all/all_rdrp_contigs.fasta",
    ("all",                "combined_full"): "result/03_RDRP_identification/10_final/all/combined_full.faa",
    ("nr_filtered",        "rdrp_contigs"): "result/03_RDRP_identification/10_final/nr_filtered/nr_filtered_contigs.fasta",
    ("nr_filtered",        "combined_full"): "result/03_RDRP_identification/10_final/nr_filtered/combined_full.faa",
    ("nr_filtered_cluster", "rdrp_contigs"): "result/03_RDRP_identification/10_final/nr_filtered_cluster/nr_filtered_cluster_contigs.fasta",
    ("nr_filtered_cluster", "combined_full"): "result/03_RDRP_identification/10_final/nr_filtered_cluster/combined_full_c90.faa",
}
_RDRP_STATS_OUTDIR = "result/99_tables/01_cp_table/03_RDRP_identification/10_final"

# ── Merged protein/contig/cluster summary tables ──────────────────────────────
_MERGED_OUTDIR         = "result/99_tables/3_RDRP_Summary"
_MERGED_PROTEIN_CONTIG = _MERGED_OUTDIR + "/nr_filtered_protein_contig_cluster.tsv"


rule all:
    input:
        expand(_OUTDIR + "/{name}", name=_TABLES.keys()),
        _ESVIRTU_INFO_TABLE,
        _RNA_VIRUS_TABLE,
        _RECALC_COUNT,
        _RECALC_RPKMF,
        _RECALC_TPM,
        _TAX_COUNT,
        _TAX_RPKMF,
        _TAX_TPM,
        [
            _RDRP_STATS_OUTDIR + f"/{stage}/{seqset}.fx2tab.tsv"
            for stage, seqset in _RDRP_STAGE_FASTA
        ],
        [
            _RDRP_STATS_OUTDIR + f"/{stage}/{seqset}.stats.tsv"
            for stage, seqset in _RDRP_STAGE_FASTA
        ],
        _MERGED_PROTEIN_CONTIG,


rule copy_table:
    input:
        lambda wildcards: ancient(_TABLES[wildcards.name])
    output:
        _OUTDIR + "/{name}"
    log:
        err = _LOG_DIR + "/{name}.err"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        "mkdir -p $(dirname {output}) && cp {input} {output} 2> {log.err}"


# ── ESvirtu info table ────────────────────────────────────────────────────────

rule esvirtu_info_table:
    input:
        rpkmf     = ancient(_OUTDIR + "/04_modified_esvirtue_ribodetector/esvirtu_rpkmf.subspecies.tsv"),
        count_tsv = ancient(_OUTDIR + "/04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv"),
        all_meta  = ancient(config["virus_pathogen_db"]["all_metadata"]),
        vhd       = ancient(config["virushostdb"]["tsv"]),
        host_meta = ancient(config["virus_pathogen_db"]["tsv"])
    output:
        _ESVIRTU_INFO_TABLE
    log:
        log = _LOG_DIR + "/esvirtu_info_table.log",
        err = _LOG_DIR + "/esvirtu_info_table.err"
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
        "python scripts/99_tables/esvirtu_info_table.py > {log.log} 2> {log.err}"


# ── RNA virus contig table (bait recovery dropped; contig_tag = RdRp evidence
# from final_rdrp_merged.tsv first, GeNomad genome-composition fallback for
# contigs with no RdRp hit) ────────────────────────────────────────────────────

rule rna_virus_contig_table:
    input:
        count_tsv   = ancient(_OUTDIR + "/04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv"),
        rdrp_merged = ancient(_OUTDIR + "/03_RDRP_identification/10_final/final_rdrp_merged.tsv"),
        genomad     = ancient(_OUTDIR + "/08_other_virus_identification/all_samples_virus_filter_summary.annotated.tsv"),
    output:
        _RNA_VIRUS_TABLE
    log:
        log = _LOG_DIR + "/rna_virus_contig_table.log",
        err = _LOG_DIR + "/rna_virus_contig_table.err"
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
        "python scripts/99_tables/rna_virus_contig_table.py > {log.log} 2> {log.err}"


# ── Recalculated RNA-virus-only abundance (esvirtu ssRNA/dsRNA rows +
# RdRpCATCH/LucaProt contig rows, RPKMF and TPM recomputed from combined
# read counts) ─────────────────────────────────────────────────────────────────

rule rna_virus_recalc_rpkmf:
    input:
        esvirtu = _ESVIRTU_INFO_TABLE,
        contig  = _RNA_VIRUS_TABLE,
    output:
        count_tsv = _RECALC_COUNT,
        rpkmf     = _RECALC_RPKMF,
        tpm       = _RECALC_TPM,
    log:
        log = _LOG_DIR + "/rna_virus_recalc_rpkmf.log",
        err = _LOG_DIR + "/rna_virus_recalc_rpkmf.err"
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
        "python scripts/99_tables/rna_virus_recalc_rpkmf.py > {log.log} 2> {log.err}"


# ── seqkit fx2tab (per-sequence: name, length, GC%) + stats (per-file summary)
# for the 03_RDRP_identification 10_final all/nr_filtered/nr_filtered_cluster
# contig + protein FASTA outputs ──────────────────────────────────────────────

rule rdrp_fasta_fx2tab:
    input:
        fasta = lambda wc: ancient(_RDRP_STAGE_FASTA[(wc.stage, wc.seqset)]),
    output:
        tsv = _RDRP_STATS_OUTDIR + "/{stage}/{seqset}.fx2tab.tsv",
    wildcard_constraints:
        stage  = "all|nr_filtered|nr_filtered_cluster",
        seqset = "rdrp_contigs|combined_full",
    log:
        err = "log/99_tables/03_RDRP_identification/10_final/{stage}/{seqset}.fx2tab.err",
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        mem_mb_per_cpu  = config["seqkit"]["memory"],
        runtime         = config["seqkit"]["runtime"],
        cpus_per_task   = config["seqkit"]["threads"],
        slurm_partition = config["seqkit"]["partition"],
        slurm_account   = config["seqkit"]["account"],
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log.err})
        echo -e "seq_id\tlength\tGC" > {output.tsv}
        seqkit fx2tab \
            --name --only-id --length --gc \
            --threads {threads} \
            {input.fasta} \
            >> {output.tsv} \
            2> {log.err}
        """


rule rdrp_fasta_stats:
    input:
        fasta = lambda wc: ancient(_RDRP_STAGE_FASTA[(wc.stage, wc.seqset)]),
    output:
        tsv = _RDRP_STATS_OUTDIR + "/{stage}/{seqset}.stats.tsv",
    wildcard_constraints:
        stage  = "all|nr_filtered|nr_filtered_cluster",
        seqset = "rdrp_contigs|combined_full",
    log:
        err = "log/99_tables/03_RDRP_identification/10_final/{stage}/{seqset}.stats.err",
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        mem_mb_per_cpu  = config["seqkit"]["memory"],
        runtime         = config["seqkit"]["runtime"],
        cpus_per_task   = config["seqkit"]["threads"],
        slurm_partition = config["seqkit"]["partition"],
        slurm_account   = config["seqkit"]["account"],
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log.err})
        seqkit stats -a --tabular \
            --threads {threads} \
            {input.fasta} \
            > {output.tsv} \
            2> {log.err}
        """


# ── Merge protein/contig/cluster info into one summary table ─────────────────
# Uses the 01_cp_table copies as input (not the original 03_RDRP_identification
# result paths) so this rule only depends on files already staged into
# 99_tables, per the user's instruction.

rule merge_rdrp_cluster_info:
    input:
        protein     = ancient(_RDRP_STATS_OUTDIR + "/nr_filtered/combined_full.fx2tab.tsv"),
        contig      = ancient(_RDRP_STATS_OUTDIR + "/nr_filtered/rdrp_contigs.fx2tab.tsv"),
        clstr       = ancient(_OUTDIR + "/03_RDRP_identification/10_final/nr_filtered_cluster/combined_full_c90.faa.clstr"),
        tophit      = ancient(_OUTDIR + "/03_RDRP_identification/21_taxonomic_assignment_cluster/full_length/full_length_diamond_rvmt_tophit.tsv"),
        rdrp_merged = ancient(_OUTDIR + "/03_RDRP_identification/10_final/final_rdrp_merged.tsv"),
    output:
        protein_contig = _MERGED_PROTEIN_CONTIG,
    log:
        log = "log/99_tables/3_RDRP_Summary/merge_rdrp_cluster_info.log",
        err = "log/99_tables/3_RDRP_Summary/merge_rdrp_cluster_info.err",
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
        "python scripts/99_tables/merge_rdrp_cluster_info.py > {log.log} 2> {log.err}"


# ── Annotate recalculated RNA-virus tables with RVMT taxonomy (contig rows
# only; esvirtu rows keep their own existing taxonomy untouched) ─────────────

rule annotate_rna_virus_taxonomy:
    input:
        count_tsv      = _RECALC_COUNT,
        rpkmf          = _RECALC_RPKMF,
        tpm            = _RECALC_TPM,
        protein_contig = ancient(_MERGED_PROTEIN_CONTIG),
    output:
        count_tsv = _TAX_COUNT,
        rpkmf     = _TAX_RPKMF,
        tpm       = _TAX_TPM,
    log:
        log = _LOG_DIR + "/annotate_rna_virus_taxonomy.log",
        err = _LOG_DIR + "/annotate_rna_virus_taxonomy.err"
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
        "python scripts/99_tables/annotate_rna_virus_taxonomy.py > {log.log} 2> {log.err}"
