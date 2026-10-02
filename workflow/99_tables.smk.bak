# Summary table generation workflow
# Collects key outputs from upstream workflows and produces final publication tables.
# Run after all upstream workflows are complete.

configfile: "config/config.yaml"

# ── Input → output mapping ────────────────────────────────────────────────────

_TABLES = {
    # RdRp identification
    "03_RDRP_identification/all_samples_rdrp_merged.tsv":
        "result/03_RDRP_identification/10_final/final_rdrp_merged.tsv",
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
    # Bait recovery
    "06_contig_bait/bait_recovered.txt":
        "result/06_contig_bait/5_bait_contigs/bait_recovered.txt",
    "06_contig_bait/bait_target_map.tsv":
        "result/06_contig_bait/4_filtered_hits/bait_target_map.tsv",
}

_OUTDIR  = "result/99_tables/01_cp_table"
_LOG_DIR = "log/99_tables"

_ESVIRTU_INFO_TABLE = "result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv"
_RNA_VIRUS_TABLE    = "result/99_tables/3_RNA_virus_table/rna_virus_contig_table.tsv"
_RECALC_RPKMF       = "result/99_tables/5_Recalculated_RPKM_result/rna_virus_recalc_rpkmf.tsv"
_RECALC_COUNT       = "result/99_tables/5_Recalculated_RPKM_result/rna_virus_recalc_read_count.tsv"
_RECALC_TPM         = "result/99_tables/5_Recalculated_RPKM_result/rna_virus_recalc_tpm.tsv"

_CORR_OUTDIR     = "result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation"
_CORR_OUTDIR_TPM = "result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation_TPM"
_CORR_THRESHOLDS = ["0.6", "0.7", "0.8", "0.9"]
_CORR_PAIR_TYPES = ["known_known_pair", "known_unknown_pair"]

_SPARCC_OUTDIR = "result/99_tables/6_Recalculated_correlation/3_sparcc"
_RUN_SPARCC   = config.get("fastspar", {}).get("run", False)

_READ_STATS_OUTDIR  = "result/99_tables/7_read_stats"
_READ_STATS_TABLE   = _READ_STATS_OUTDIR + "/sample_read_stats.tsv"
_FASTP_DIR          = config["rna_fastp_dir"]
_RIBODETECTOR_DIR   = config["ribodetector_dir"]
_SAMPLES            = config["rna_samples"]


rule all:
    input:
        expand(_OUTDIR + "/{name}", name=_TABLES.keys()),
        _ESVIRTU_INFO_TABLE,
        _RNA_VIRUS_TABLE,
        _RECALC_RPKMF,
        _RECALC_COUNT,
        _RECALC_TPM,
        expand(
            _CORR_OUTDIR + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ),
        expand(
            _CORR_OUTDIR + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
            threshold=_CORR_THRESHOLDS,
        ),
        expand(
            _CORR_OUTDIR_TPM + "/2_spearman_analysis/{pair_type}/r{threshold}.tsv",
            pair_type=_CORR_PAIR_TYPES,
            threshold=_CORR_THRESHOLDS,
        ),
        expand(
            _CORR_OUTDIR_TPM + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
            threshold=_CORR_THRESHOLDS,
        ),
        _READ_STATS_TABLE,
        [_SPARCC_OUTDIR + "/sparcc_pairs.tsv"] if _RUN_SPARCC else [],


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


# ── RNA virus contig table ────────────────────────────────────────────────────

rule rna_virus_contig_table:
    input:
        rpkmf       = ancient(_OUTDIR + "/04_modified_esvirtue_ribodetector/esvirtu_rpkmf.subspecies.tsv"),
        count_tsv   = ancient(_OUTDIR + "/04_modified_esvirtue_ribodetector/esvirtu_read_count.subspecies.tsv"),
        rdrp_merged = ancient(_OUTDIR + "/03_RDRP_identification/all_samples_rdrp_merged.tsv"),
        bait        = ancient(_OUTDIR + "/06_contig_bait/bait_target_map.tsv"),
        genomad     = ancient(_OUTDIR + "/08_other_virus_identification/all_samples_virus_filter_summary.annotated.tsv")
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


# ── Recalculated RPKMF (RNA virus only) ──────────────────────────────────────

rule rna_virus_recalc_rpkmf:
    input:
        esvirtu = _ESVIRTU_INFO_TABLE,
        contig  = _RNA_VIRUS_TABLE
    output:
        rpkmf     = _RECALC_RPKMF,
        count_tsv = _RECALC_COUNT,
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
        all_vs_all      = _CORR_OUTDIR     + "/1_all_vs_all/rna_virus.all_vs_all.spearman.filtered.tsv",
        pair_outdir     = _CORR_OUTDIR     + "/2_spearman_analysis",
        all_vs_all_tpm  = _CORR_OUTDIR_TPM + "/1_all_vs_all/rna_virus.all_vs_all.spearman.filtered.tsv",
        pair_outdir_tpm = _CORR_OUTDIR_TPM + "/2_spearman_analysis",
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
        """


rule recalc_annotate_rdrp_category:
    input:
        ku_tsv       = _CORR_OUTDIR     + "/2_spearman_analysis/known_unknown_pair/r{threshold}.tsv",
        ku_tsv_tpm   = _CORR_OUTDIR_TPM + "/2_spearman_analysis/known_unknown_pair/r{threshold}.tsv",
        contig_table = _RNA_VIRUS_TABLE,
    output:
        annotated     = _CORR_OUTDIR     + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
        annotated_tpm = _CORR_OUTDIR_TPM + "/3_spearman_analysis_RdRp_annotated/known_unknown_pair/r{threshold}.annotated.tsv",
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
        """


# ── Per-sample read statistics (after fastp and after ribodetector) ───────────

rule sample_read_stats:
    input:
        fastp_r1 = expand(_FASTP_DIR + "/{sample}/{sample}_1P.fq.gz", sample=_SAMPLES),
        fastp_r2 = expand(_FASTP_DIR + "/{sample}/{sample}_2P.fq.gz", sample=_SAMPLES),
        ribo_r1  = expand(_RIBODETECTOR_DIR + "/{sample}/{sample}_nonrrna.1.fq.gz", sample=_SAMPLES),
        ribo_r2  = expand(_RIBODETECTOR_DIR + "/{sample}/{sample}_nonrrna.2.fq.gz", sample=_SAMPLES),
    output:
        _READ_STATS_TABLE,
    log:
        out = _LOG_DIR + "/sample_read_stats.log",
        err = _LOG_DIR + "/sample_read_stats.err",
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    params:
        fastp_dir = _FASTP_DIR,
        ribo_dir  = _RIBODETECTOR_DIR,
        samples   = _SAMPLES,
    shell:
        """
        mkdir -p {_READ_STATS_OUTDIR} log/99_tables

        # header
        echo -e "sample\tfastp_R1\tfastp_R2\tfastp_total\tribodetector_R1\tribodetector_R2\tribodetector_total" > {output}

        for sample in {params.samples}; do
            fastp_r1={params.fastp_dir}/${{sample}}/${{sample}}_1P.fq.gz
            fastp_r2={params.fastp_dir}/${{sample}}/${{sample}}_2P.fq.gz
            ribo_r1={params.ribo_dir}/${{sample}}/${{sample}}_nonrrna.1.fq.gz
            ribo_r2={params.ribo_dir}/${{sample}}/${{sample}}_nonrrna.2.fq.gz

            fp_r1=$(seqkit stats -j {threads} "$fastp_r1" 2>/dev/null | awk 'NR==2 {{print $4}}' | tr -d ',')
            fp_r2=$(seqkit stats -j {threads} "$fastp_r2" 2>/dev/null | awk 'NR==2 {{print $4}}' | tr -d ',')
            fp_total=$((fp_r1 + fp_r2))

            rb_r1=$(seqkit stats -j {threads} "$ribo_r1" 2>/dev/null | awk 'NR==2 {{print $4}}' | tr -d ',')
            rb_r2=$(seqkit stats -j {threads} "$ribo_r2" 2>/dev/null | awk 'NR==2 {{print $4}}' | tr -d ',')
            rb_total=$((rb_r1 + rb_r2))

            echo -e "${{sample}}\t${{fp_r1}}\t${{fp_r2}}\t${{fp_total}}\t${{rb_r1}}\t${{rb_r2}}\t${{rb_total}}"
        done >> {output}

        echo "Done: $(wc -l < {output}) samples" > {log.out}
        """


# ── SparCC (via FastSpar) correlation on read counts ─────────────────────────

if _RUN_SPARCC:

    rule sparcc_prepare:
        input:
            counts = ancient(_RECALC_COUNT),
            script = config["scripts"]["prepare_fastspar_input"],
        output:
            otu = _SPARCC_OUTDIR + "/otu_table.tsv",
        params:
            min_samples_others = config["fastspar"]["min_samples_others"],
        log:
            err = _LOG_DIR + "/sparcc/prepare.err",
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
            mkdir -p {_SPARCC_OUTDIR} log/99_tables/sparcc
            python {input.script} \
                --input              {input.counts} \
                --output             {output.otu} \
                --min-samples-others {params.min_samples_others} \
                2> {log.err}
            """

    rule sparcc_run:
        input:
            otu = _SPARCC_OUTDIR + "/otu_table.tsv",
        output:
            corr = _SPARCC_OUTDIR + "/correlation.tsv",
            cov  = _SPARCC_OUTDIR + "/covariance.tsv",
            pval = _SPARCC_OUTDIR + "/pvalues.tsv",
        params:
            outdir     = _SPARCC_OUTDIR,
            iterations = config["fastspar"]["iterations"],
            exclude_it = config["fastspar"]["exclude_iterations"],
            threshold  = config["fastspar"]["threshold"],
            bootstraps = config["fastspar"]["bootstraps"],
        log:
            out = _LOG_DIR + "/sparcc/fastspar.log",
            err = _LOG_DIR + "/sparcc/fastspar.err",
        conda:
            "../envs/fastspar.yaml"
        threads: config["fastspar"]["threads"]
        resources:
            mem_mb_per_cpu  = config["fastspar"]["memory"],
            runtime         = config["fastspar"]["runtime"],
            cpus_per_task   = config["fastspar"]["threads"],
            slurm_partition = config["fastspar"]["partition"],
            slurm_account   = config["fastspar"]["account"],
        shell:
            """
            mkdir -p {params.outdir}/bootstraps_counts \
                     {params.outdir}/bootstraps_correlation \
                     log/99_tables/sparcc

            fastspar \
                --otu_table  {input.otu} \
                --correlation {output.corr} \
                --covariance  {output.cov} \
                --iterations  {params.iterations} \
                -x {params.exclude_it} \
                --threshold   {params.threshold} \
                --threads     {threads} \
                > {log.out} 2> {log.err}

            fastspar_bootstrap \
                --otu_table  {input.otu} \
                --number     {params.bootstraps} \
                --prefix     {params.outdir}/bootstraps_counts/bootstrap \
                --threads    {threads} \
                >> {log.out} 2>> {log.err}

            # Bootstrap correlations: sequential, each fastspar uses --threads internally
            for f in {params.outdir}/bootstraps_counts/bootstrap_*.tsv; do
                base=$(basename "$f" .tsv)
                fastspar \
                    --otu_table "$f" \
                    --correlation {params.outdir}/bootstraps_correlation/cor_${{base}}.tsv \
                    --covariance  {params.outdir}/bootstraps_correlation/cov_${{base}}.tsv \
                    --iterations {params.iterations} \
                    -x {params.exclude_it} \
                    --threshold {params.threshold} \
                    --threads {threads} --yes \
                    >> {log.out} 2>> {log.err}
            done

            fastspar_pvalues \
                --otu_table      {input.otu} \
                --correlation    {output.corr} \
                --prefix         {params.outdir}/bootstraps_correlation/cor_bootstrap_ \
                --permutations   {params.bootstraps} \
                --outfile        {output.pval} \
                --threads        {threads} \
                >> {log.out} 2>> {log.err}
            """

    rule sparcc_parse:
        input:
            corr   = _SPARCC_OUTDIR + "/correlation.tsv",
            pval   = _SPARCC_OUTDIR + "/pvalues.tsv",
            script = config["scripts"]["parse_fastspar_output"],
        output:
            pairs = _SPARCC_OUTDIR + "/sparcc_pairs.tsv",
        params:
            thresholds  = config["fastspar"]["r_threshold"],
            r_min       = config["fastspar"]["r_min"],
            p_threshold = config["fastspar"]["p_threshold"],
        log:
            err = _LOG_DIR + "/sparcc/parse.err",
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
            python {input.script} \
                --correlation {input.corr} \
                --pvalues     {input.pval} \
                --output      {output.pairs} \
                --thresholds  {params.thresholds} \
                --r-threshold {params.r_min} \
                --p-threshold {params.p_threshold} \
                2> {log.err}
            """
