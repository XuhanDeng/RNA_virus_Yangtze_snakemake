# Summary table generation workflow
# Collects key outputs from upstream workflows and produces final publication tables.
# Run after all upstream workflows are complete.

configfile: "config/config.yaml"

# ── Input → output mapping ────────────────────────────────────────────────────

_TABLES = {
    # RdRp identification
    "03_RDRP_identification/tier_summary.tsv":
        "result/03_RDRP_identification/4_motif_search/motif_results/tier_summary.tsv",
    "03_RDRP_identification/motif_hits_best.tsv":
        "result/03_RDRP_identification/4_motif_search/motif_results/motif_hits_best.tsv",
    "03_RDRP_identification/motif_order_complete.tsv":
        "result/03_RDRP_identification/4_motif_search/motif_results/motif_order_complete.tsv",
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
    "04_03_esvirtu_rdrp_identification/esvirtu_reference_metadata.tsv":
        "result/04_03_esvirtu_rdrp_identification/2_metadata/esvirtu_reference.metadata.tsv",
    # Other virus identification (GeNomad + ICTV VMR annotation)
    "08_other_virus_identification/all_samples_virus_filter_summary.annotated.tsv":
        "result/08_other_virus_identification/2_virus_filter_summary/all_samples_virus_filter_summary.annotated.tsv",
    # Bait recovery
    "06_contig_bait/bait_recovered.txt":
        "result/06_contig_bait/5_bait_contigs/bait_recovered.txt",
    "06_contig_bait/bait_target_map.tsv":
        "result/06_contig_bait/4_filtered_hits/bait_target_map.tsv",
    # Correlation -- known-known pairs (subspecies)
    "05_esvirtu_correlation_analysis/known_known_pair/r0.6.tsv":
        "result/05_esvirtu_correlation_analysis/2_spearman_analysis/subspecies/known_known_pair/r0.6.tsv",
    "05_esvirtu_correlation_analysis/known_known_pair/r0.7.tsv":
        "result/05_esvirtu_correlation_analysis/2_spearman_analysis/subspecies/known_known_pair/r0.7.tsv",
    "05_esvirtu_correlation_analysis/known_known_pair/r0.8.tsv":
        "result/05_esvirtu_correlation_analysis/2_spearman_analysis/subspecies/known_known_pair/r0.8.tsv",
    "05_esvirtu_correlation_analysis/known_known_pair/r0.9.tsv":
        "result/05_esvirtu_correlation_analysis/2_spearman_analysis/subspecies/known_known_pair/r0.9.tsv",
    # Correlation -- known-unknown pairs annotated (subspecies)
    "05_esvirtu_correlation_analysis/known_unknown_pair/r0.6.annotated.tsv":
        "result/05_esvirtu_correlation_analysis/3_spearman_analysis_RdRp_annotated/subspecies/known_unknown_pair/r0.6.annotated.tsv",
    "05_esvirtu_correlation_analysis/known_unknown_pair/r0.7.annotated.tsv":
        "result/05_esvirtu_correlation_analysis/3_spearman_analysis_RdRp_annotated/subspecies/known_unknown_pair/r0.7.annotated.tsv",
    "05_esvirtu_correlation_analysis/known_unknown_pair/r0.8.annotated.tsv":
        "result/05_esvirtu_correlation_analysis/3_spearman_analysis_RdRp_annotated/subspecies/known_unknown_pair/r0.8.annotated.tsv",
    "05_esvirtu_correlation_analysis/known_unknown_pair/r0.9.annotated.tsv":
        "result/05_esvirtu_correlation_analysis/3_spearman_analysis_RdRp_annotated/subspecies/known_unknown_pair/r0.9.annotated.tsv",
}

_OUTDIR  = "result/99_tables/01_cp_table"
_LOG_DIR = "log/99_tables"


_ESVIRTU_INFO_TABLE  = "result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv"
_RNA_VIRUS_TABLE     = "result/99_tables/3_RNA_virus_table/rna_virus_contig_table.tsv"
_RECALC_RPKMF        = "result/99_tables/5_Recalculated_RPKM_result/1_RNA_virus_Recalculated_RPKM/rna_virus_recalc_rpkmf.tsv"
_RECALC_COUNT        = "result/99_tables/5_Recalculated_RPKM_result/1_RNA_virus_Recalculated_RPKM/rna_virus_recalc_read_count.tsv"
_RECALC_TPM          = "result/99_tables/5_Recalculated_RPKM_result/1_RNA_virus_Recalculated_RPKM/rna_virus_recalc_tpm.tsv"

_CORR_SCRIPTS    = "scripts/99_tables/06_recalc_correlation"
_CORR_OUTDIR     = "result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation"
_CORR_OUTDIR_TPM = "result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation_TPM"
_CORR_THRESHOLDS = ["0.6", "0.7", "0.8", "0.9"]
_CORR_PAIR_TYPES = ["known_known_pair", "known_unknown_pair"]

_RECALC_RPKMF_R  = "result/99_tables/5_Recalculated_RPKM_result/2_RNA_virus_Recalculated_RPKM_restricted/rna_virus_recalc_rpkmf.tsv"
_RECALC_COUNT_R  = "result/99_tables/5_Recalculated_RPKM_result/2_RNA_virus_Recalculated_RPKM_restricted/rna_virus_recalc_read_count.tsv"
_RECALC_TPM_R    = "result/99_tables/5_Recalculated_RPKM_result/2_RNA_virus_Recalculated_RPKM_restricted/rna_virus_recalc_tpm.tsv"
_CORR_OUTDIR_R   = "result/99_tables/6_Recalculated_correlation/2_RNA_virus_Recalculated_correlation_restricted"
_CORR_OUTDIR_R_TPM = "result/99_tables/6_Recalculated_correlation/2_RNA_virus_Recalculated_correlation_restricted_TPM"

_RECALC_RPKMF_R2 = "result/99_tables/5_Recalculated_RPKM_result/3_RNA_virus_Recalculated_RPKM_tier1only/rna_virus_recalc_rpkmf.tsv"
_RECALC_COUNT_R2 = "result/99_tables/5_Recalculated_RPKM_result/3_RNA_virus_Recalculated_RPKM_tier1only/rna_virus_recalc_read_count.tsv"
_RECALC_TPM_R2   = "result/99_tables/5_Recalculated_RPKM_result/3_RNA_virus_Recalculated_RPKM_tier1only/rna_virus_recalc_tpm.tsv"
_CORR_OUTDIR_R2  = "result/99_tables/6_Recalculated_correlation/3_RNA_virus_Recalculated_correlation_tier1only"
_CORR_OUTDIR_R2_TPM = "result/99_tables/6_Recalculated_correlation/3_RNA_virus_Recalculated_correlation_tier1only_TPM"

_RUN_RESTRICTED  = config.get("recalc_rpkmf", {}).get("run_restricted", False)
_RUN_TIER1_ONLY  = config.get("recalc_rpkmf", {}).get("run_tier1_only", False)

rule all:
    input:
        expand(_OUTDIR + "/{name}", name=_TABLES.keys()),
        _ESVIRTU_INFO_TABLE,
        _RNA_VIRUS_TABLE,
        _RECALC_RPKMF,
        _RECALC_COUNT,
        expand(
            _CORR_OUTDIR + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ),
        expand(
            _CORR_OUTDIR + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
            threshold=_CORR_THRESHOLDS,
        ),
        _RECALC_TPM,
        expand(
            _CORR_OUTDIR_TPM + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ),
        expand(
            _CORR_OUTDIR_TPM + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
            threshold=_CORR_THRESHOLDS,
        ),
        *([
            _RECALC_RPKMF_R,
            _RECALC_COUNT_R,
            _RECALC_TPM_R,
            expand(
                _CORR_OUTDIR_R + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
                pair_type=_CORR_PAIR_TYPES,
                threshold=_CORR_THRESHOLDS,
            ),
            expand(
                _CORR_OUTDIR_R + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
                threshold=_CORR_THRESHOLDS,
            ),
            expand(
                _CORR_OUTDIR_R_TPM + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
                pair_type=_CORR_PAIR_TYPES,
                threshold=_CORR_THRESHOLDS,
            ),
            expand(
                _CORR_OUTDIR_R_TPM + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
                threshold=_CORR_THRESHOLDS,
            ),
        ] if _RUN_RESTRICTED else []),
        *([
            _RECALC_RPKMF_R2,
            _RECALC_COUNT_R2,
            _RECALC_TPM_R2,
            expand(
                _CORR_OUTDIR_R2 + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
                pair_type=_CORR_PAIR_TYPES,
                threshold=_CORR_THRESHOLDS,
            ),
            expand(
                _CORR_OUTDIR_R2 + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
                threshold=_CORR_THRESHOLDS,
            ),
            expand(
                _CORR_OUTDIR_R2_TPM + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
                pair_type=_CORR_PAIR_TYPES,
                threshold=_CORR_THRESHOLDS,
            ),
            expand(
                _CORR_OUTDIR_R2_TPM + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
                threshold=_CORR_THRESHOLDS,
            ),
        ] if _RUN_TIER1_ONLY else [])


rule copy_table:
    input:
        lambda wildcards: ancient(_TABLES[wildcards.name])
    output:
        _OUTDIR + "/{name}"
    log:
        err = _LOG_DIR + "/{name}.err"
    shell:
        "mkdir -p $(dirname {output}) && cp {input} {output} 2> {log.err}"


# ── ESvirtu info table ────────────────────────────────────────────────────────

rule esvirtu_info_table:
    input:
        rpkmf     = ancient(_TABLES["04_modified_esvirtue_ribodetector/esvirtu_rpkmf.subspecies.tsv"]),
        count_tsv = ancient(_TABLES["04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv"]),
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
    shell:
        "python scripts/99_tables/esvirtu_info_table.py > {log.log} 2> {log.err}"


# ── RNA virus contig table ────────────────────────────────────────────────────

rule rna_virus_contig_table:
    input:
        rpkmf = ancient(_OUTDIR + "/04_modified_esvirtue_ribodetector/esvirtu_rpkmf.subspecies.tsv"),
        count_tsv = ancient(_OUTDIR + "/04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv"),
        tier  = ancient(_OUTDIR + "/03_RDRP_identification/tier_summary.tsv"),
        bait  = ancient(_OUTDIR + "/06_contig_bait/bait_target_map.tsv"),
        genomad = ancient(_OUTDIR + "/08_other_virus_identification/all_samples_virus_filter_summary.annotated.tsv")
    output:
        _RNA_VIRUS_TABLE
    log:
        log = _LOG_DIR + "/rna_virus_contig_table.log",
        err = _LOG_DIR + "/rna_virus_contig_table.err"
    conda:
        "../envs/python.yaml"
    shell:
        "python scripts/99_tables/rna_virus_contig_table.py > {log.log} 2> {log.err}"


# ── Recalculated RPKMF (RNA virus only) ──────────────────────────────────────

rule rna_virus_recalc_rpkmf:
    input:
        esvirtu = _ESVIRTU_INFO_TABLE,
        contig  = _RNA_VIRUS_TABLE
    output:
        rpkmf     = _RECALC_RPKMF,
        count_tsv = _RECALC_COUNT,
        tpm       = _RECALC_TPM,
        **( {"rpkmf_r": _RECALC_RPKMF_R, "count_tsv_r": _RECALC_COUNT_R, "tpm_r": _RECALC_TPM_R}
            if _RUN_RESTRICTED else {} ),
        **( {"rpkmf_r2": _RECALC_RPKMF_R2, "count_tsv_r2": _RECALC_COUNT_R2, "tpm_r2": _RECALC_TPM_R2}
            if _RUN_TIER1_ONLY else {} )
    log:
        log = _LOG_DIR + "/rna_virus_recalc_rpkmf.log",
        err = _LOG_DIR + "/rna_virus_recalc_rpkmf.err"
    conda:
        "../envs/python.yaml"
    shell:
        """
        python scripts/99_tables/rna_virus_recalc_rpkmf.py > {log.log} 2> {log.err}
        """ + ("""
        python scripts/99_tables/rna_virus_recalc_rpkmf_restricted.py >> {log.log} 2>> {log.err}
        """ if _RUN_RESTRICTED else "") + ("""
        python scripts/99_tables/rna_virus_recalc_rpkmf_tier1only.py >> {log.log} 2>> {log.err}
        """ if _RUN_TIER1_ONLY else "")


# ── Recalculated correlation (RNA virus only) ─────────────────────────────────

rule recalc_spearman_and_filter:
    input:
        table = _RECALC_RPKMF,
        tpm   = _RECALC_TPM
    output:
        expand(
            _CORR_OUTDIR + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ),
        expand(
            _CORR_OUTDIR_TPM + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ),
        *( expand(
            _CORR_OUTDIR_R + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ) if _RUN_RESTRICTED else [] ),
        *( expand(
            _CORR_OUTDIR_R_TPM + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ) if _RUN_RESTRICTED else [] ),
        *( expand(
            _CORR_OUTDIR_R2 + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ) if _RUN_TIER1_ONLY else [] ),
        *( expand(
            _CORR_OUTDIR_R2_TPM + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ) if _RUN_TIER1_ONLY else [] )
    log:
        out = _LOG_DIR + "/6_recalc_correlation/spearman_all_vs_all.log",
        err = _LOG_DIR + "/6_recalc_correlation/spearman_all_vs_all.err"
    conda:
        "../envs/python.yaml"
    threads: config["correlation"]["threads"]
    resources:
        mem_mb_per_cpu  = config["correlation"]["memory"],
        runtime         = config["correlation"]["runtime"],
        cpus_per_task   = config["correlation"]["threads"],
        slurm_partition = config["correlation"]["partition"],
        slurm_account   = config["correlation"]["account"],
    params:
        all_vs_all  = _CORR_OUTDIR   + "/1_all_vs_all/rna_virus.all_vs_all.spearman.filtered.tsv",
        pair_outdir = _CORR_OUTDIR   + "/2_spearman_analysis",
        all_vs_all_tpm  = _CORR_OUTDIR_TPM + "/1_all_vs_all/rna_virus.all_vs_all.spearman.filtered.tsv",
        pair_outdir_tpm = _CORR_OUTDIR_TPM + "/2_spearman_analysis",
        all_vs_all_r  = _CORR_OUTDIR_R  + "/1_all_vs_all/rna_virus.all_vs_all.spearman.filtered.tsv",
        pair_outdir_r = _CORR_OUTDIR_R  + "/2_spearman_analysis",
        all_vs_all_r_tpm  = _CORR_OUTDIR_R_TPM + "/1_all_vs_all/rna_virus.all_vs_all.spearman.filtered.tsv",
        pair_outdir_r_tpm = _CORR_OUTDIR_R_TPM + "/2_spearman_analysis",
        all_vs_all_r2  = _CORR_OUTDIR_R2 + "/1_all_vs_all/rna_virus.all_vs_all.spearman.filtered.tsv",
        pair_outdir_r2 = _CORR_OUTDIR_R2 + "/2_spearman_analysis",
        all_vs_all_r2_tpm  = _CORR_OUTDIR_R2_TPM + "/1_all_vs_all/rna_virus.all_vs_all.spearman.filtered.tsv",
        pair_outdir_r2_tpm = _CORR_OUTDIR_R2_TPM + "/2_spearman_analysis"
    shell:
        """
        mkdir -p $(dirname {params.all_vs_all}) {params.pair_outdir}
        python scripts/99_tables/06_recalc_correlation/spearman_all_vs_all.py \
            --input            {input.table} \
            --output-filtered  {params.all_vs_all} \
            --min-samples      {config[correlation][min_samples]} \
            --thresholds       {config[correlation][thresholds]} \
            --p-threshold      {config[correlation][p_threshold]} \
            --chunk-size       {config[correlation][chunk_size_all_vs_all]} \
            --threads          {threads} \
            > {log.out} 2> {log.err}
        python scripts/99_tables/06_recalc_correlation/filter_spearman_pairs.py \
            --input  {params.all_vs_all} \
            --outdir {params.pair_outdir} \
            2>> {log.err}

        mkdir -p $(dirname {params.all_vs_all_tpm}) {params.pair_outdir_tpm}
        python scripts/99_tables/06_recalc_correlation/spearman_all_vs_all.py \
            --input            {input.tpm} \
            --value-suffix     _tpm \
            --output-filtered  {params.all_vs_all_tpm} \
            --min-samples      {config[correlation][min_samples]} \
            --thresholds       {config[correlation][thresholds]} \
            --p-threshold      {config[correlation][p_threshold]} \
            --chunk-size       {config[correlation][chunk_size_all_vs_all]} \
            --threads          {threads} \
            >> {log.out} 2>> {log.err}
        python scripts/99_tables/06_recalc_correlation/filter_spearman_pairs.py \
            --input  {params.all_vs_all_tpm} \
            --outdir {params.pair_outdir_tpm} \
            2>> {log.err}
        """ + ("""
        mkdir -p $(dirname {params.all_vs_all_r}) {params.pair_outdir_r}
        python scripts/99_tables/06_recalc_correlation/spearman_all_vs_all.py \
            --input            """ + _RECALC_RPKMF_R + """ \
            --output-filtered  {params.all_vs_all_r} \
            --min-samples      {config[correlation][min_samples]} \
            --thresholds       {config[correlation][thresholds]} \
            --p-threshold      {config[correlation][p_threshold]} \
            --chunk-size       {config[correlation][chunk_size_all_vs_all]} \
            --threads          {threads} \
            >> {log.out} 2>> {log.err}
        python scripts/99_tables/06_recalc_correlation/filter_spearman_pairs.py \
            --input  {params.all_vs_all_r} \
            --outdir {params.pair_outdir_r} \
            2>> {log.err}

        mkdir -p $(dirname {params.all_vs_all_r_tpm}) {params.pair_outdir_r_tpm}
        python scripts/99_tables/06_recalc_correlation/spearman_all_vs_all.py \
            --input            """ + _RECALC_TPM_R + """ \
            --value-suffix     _tpm \
            --output-filtered  {params.all_vs_all_r_tpm} \
            --min-samples      {config[correlation][min_samples]} \
            --thresholds       {config[correlation][thresholds]} \
            --p-threshold      {config[correlation][p_threshold]} \
            --chunk-size       {config[correlation][chunk_size_all_vs_all]} \
            --threads          {threads} \
            >> {log.out} 2>> {log.err}
        python scripts/99_tables/06_recalc_correlation/filter_spearman_pairs.py \
            --input  {params.all_vs_all_r_tpm} \
            --outdir {params.pair_outdir_r_tpm} \
            2>> {log.err}
        """ if _RUN_RESTRICTED else "") + ("""
        mkdir -p $(dirname {params.all_vs_all_r2}) {params.pair_outdir_r2}
        python scripts/99_tables/06_recalc_correlation/spearman_all_vs_all.py \
            --input            """ + _RECALC_RPKMF_R2 + """ \
            --output-filtered  {params.all_vs_all_r2} \
            --min-samples      {config[correlation][min_samples]} \
            --thresholds       {config[correlation][thresholds]} \
            --p-threshold      {config[correlation][p_threshold]} \
            --chunk-size       {config[correlation][chunk_size_all_vs_all]} \
            --threads          {threads} \
            >> {log.out} 2>> {log.err}
        python scripts/99_tables/06_recalc_correlation/filter_spearman_pairs.py \
            --input  {params.all_vs_all_r2} \
            --outdir {params.pair_outdir_r2} \
            2>> {log.err}

        mkdir -p $(dirname {params.all_vs_all_r2_tpm}) {params.pair_outdir_r2_tpm}
        python scripts/99_tables/06_recalc_correlation/spearman_all_vs_all.py \
            --input            """ + _RECALC_TPM_R2 + """ \
            --value-suffix     _tpm \
            --output-filtered  {params.all_vs_all_r2_tpm} \
            --min-samples      {config[correlation][min_samples]} \
            --thresholds       {config[correlation][thresholds]} \
            --p-threshold      {config[correlation][p_threshold]} \
            --chunk-size       {config[correlation][chunk_size_all_vs_all]} \
            --threads          {threads} \
            >> {log.out} 2>> {log.err}
        python scripts/99_tables/06_recalc_correlation/filter_spearman_pairs.py \
            --input  {params.all_vs_all_r2_tpm} \
            --outdir {params.pair_outdir_r2_tpm} \
            2>> {log.err}
        """ if _RUN_TIER1_ONLY else "")


rule recalc_annotate_rdrp_category:
    input:
        ku_tsv        = _CORR_OUTDIR + "/2_spearman_analysis/known_unknown_pair/r{threshold}.tsv",
        ku_tsv_tpm    = _CORR_OUTDIR_TPM + "/2_spearman_analysis/known_unknown_pair/r{threshold}.tsv",
        contig_table  = _RNA_VIRUS_TABLE,
        **( {"ku_tsv_r": _CORR_OUTDIR_R + "/2_spearman_analysis/known_unknown_pair/r{threshold}.tsv"}
            if _RUN_RESTRICTED else {} ),
        **( {"ku_tsv_r_tpm": _CORR_OUTDIR_R_TPM + "/2_spearman_analysis/known_unknown_pair/r{threshold}.tsv"}
            if _RUN_RESTRICTED else {} ),
        **( {"ku_tsv_r2": _CORR_OUTDIR_R2 + "/2_spearman_analysis/known_unknown_pair/r{threshold}.tsv"}
            if _RUN_TIER1_ONLY else {} ),
        **( {"ku_tsv_r2_tpm": _CORR_OUTDIR_R2_TPM + "/2_spearman_analysis/known_unknown_pair/r{threshold}.tsv"}
            if _RUN_TIER1_ONLY else {} )
    output:
        annotated = _CORR_OUTDIR + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
        annotated_tpm = _CORR_OUTDIR_TPM + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
        **( {"annotated_r": _CORR_OUTDIR_R + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv"}
            if _RUN_RESTRICTED else {} ),
        **( {"annotated_r_tpm": _CORR_OUTDIR_R_TPM + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv"}
            if _RUN_RESTRICTED else {} ),
        **( {"annotated_r2": _CORR_OUTDIR_R2 + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv"}
            if _RUN_TIER1_ONLY else {} ),
        **( {"annotated_r2_tpm": _CORR_OUTDIR_R2_TPM + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv"}
            if _RUN_TIER1_ONLY else {} )
    wildcard_constraints:
        threshold = "0\\.6|0\\.7|0\\.8|0\\.9"
    log:
        err = _LOG_DIR + "/6_recalc_correlation/annotate_r{threshold}.err"
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
        """
        mkdir -p $(dirname {output.annotated})
        python scripts/99_tables/06_recalc_correlation/annotate_rdrp_category.py \
            --input          {input.ku_tsv} \
            --output         {output.annotated} \
            --contig-table   {input.contig_table} \
            2> {log.err}

        mkdir -p $(dirname {output.annotated_tpm})
        python scripts/99_tables/06_recalc_correlation/annotate_rdrp_category.py \
            --input          {input.ku_tsv_tpm} \
            --output         {output.annotated_tpm} \
            --contig-table   {input.contig_table} \
            2>> {log.err}
        """ + ("""
        mkdir -p $(dirname {output.annotated_r})
        python scripts/99_tables/06_recalc_correlation/annotate_rdrp_category.py \
            --input          {input.ku_tsv_r} \
            --output         {output.annotated_r} \
            --contig-table   {input.contig_table} \
            2>> {log.err}

        mkdir -p $(dirname {output.annotated_r_tpm})
        python scripts/99_tables/06_recalc_correlation/annotate_rdrp_category.py \
            --input          {input.ku_tsv_r_tpm} \
            --output         {output.annotated_r_tpm} \
            --contig-table   {input.contig_table} \
            2>> {log.err}
        """ if _RUN_RESTRICTED else "") + ("""
        mkdir -p $(dirname {output.annotated_r2})
        python scripts/99_tables/06_recalc_correlation/annotate_rdrp_category.py \
            --input          {input.ku_tsv_r2} \
            --output         {output.annotated_r2} \
            --contig-table   {input.contig_table} \
            2>> {log.err}

        mkdir -p $(dirname {output.annotated_r2_tpm})
        python scripts/99_tables/06_recalc_correlation/annotate_rdrp_category.py \
            --input          {input.ku_tsv_r2_tpm} \
            --output         {output.annotated_r2_tpm} \
            --contig-table   {input.contig_table} \
            2>> {log.err}
        """ if _RUN_TIER1_ONLY else "")
