configfile: "config/config.yaml"


_UPSTREAM = "result/04_modified_esvirtue_ribodetector/es/Merge"

# input read_count tables (raw counts — required by FastSpar)
_RC_INPUT = {
    "assembly":   f"{_UPSTREAM}/all_samples.detected_virus.assembly_summary.read_count.tsv",
    "subspecies": f"{_UPSTREAM}/all_samples.detected_virus.assembly_summary.subspecies.read_count.tsv",
    "species":    f"{_UPSTREAM}/all_samples.detected_virus.assembly_summary.species.read_count.tsv",
}

_LEVELS = list(_RC_INPUT.keys())

# _OUTDIR_SPARCC  = "result/05_esvirtu_correlation_analysis/3_fastspar"   # commented out
_OUTDIR_ALL2ALL   = "result/05_esvirtu_correlation_analysis/2_all_vs_all"
_OUTDIR_ANALYSIS  = "result/05_esvirtu_correlation_analysis/3_spearman_analysis"

_FILTER_LEVELS    = ["subspecies", "species"]
_THRESHOLDS       = ["0.6", "0.7", "0.8", "0.9"]
_PAIR_TYPES       = ["known_known_pair", "known_unknown_pair"]

_TIER_SUMMARY_TSV     = "result/03_RDRP_identification/4_motif_search/motif_results/tier_summary.tsv"

_OUTDIR_ANNOTATED = "result/05_esvirtu_correlation_analysis/4_spearman_analysis_RdRp_annotated"


rule all:
    input:
        expand(f"{_OUTDIR_ALL2ALL}/{{level}}.all_vs_all.spearman.filtered.tsv", level=_LEVELS),
        # FastSpar outputs — commented out
        # expand(f"{_OUTDIR_SPARCC}/{{level}}/{{level}}.sparcc.filtered.tsv", level=_LEVELS),
        expand(
            f"{_OUTDIR_ANALYSIS}/{{level}}/{{pair_type}}/r{{threshold}}.tsv",
            level=_FILTER_LEVELS,
            pair_type=_PAIR_TYPES,
            threshold=_THRESHOLDS,
        ),
        expand(
            f"{_OUTDIR_ANNOTATED}/{{level}}/known_unknown_pair/r{{threshold}}.annotated.tsv",
            level=_FILTER_LEVELS,
            threshold=_THRESHOLDS,
        ),


# ------------------------------------------------------------------ #
# Spearman all-vs-all (raw RPKMF, existing script)                   #
# ------------------------------------------------------------------ #

_RPKMF_INPUT = {
    "assembly":   f"{_UPSTREAM}/all_samples.detected_virus.assembly_summary.rpkmf.tsv",
    "subspecies": f"{_UPSTREAM}/all_samples.detected_virus.assembly_summary.rpkmf.subspecies.tsv",
    "species":    f"{_UPSTREAM}/all_samples.detected_virus.assembly_summary.rpkmf.species.tsv",
}

rule spearman_all_vs_all:
    input:
        table = lambda wc: ancient(_RPKMF_INPUT[wc.level]),
    output:
        filt_corr = f"{_OUTDIR_ALL2ALL}/{{level}}.all_vs_all.spearman.filtered.tsv",
    wildcard_constraints:
        level = "assembly|subspecies|species",
    threads: config["correlation"]["threads"]
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb_per_cpu  = config["correlation"]["memory"],
        runtime         = config["correlation"]["runtime"],
        cpus_per_task   = config["correlation"]["threads"],
        slurm_partition = config["correlation"]["partition"],
        slurm_account   = config["correlation"]["account"],
    params:
        script      = config["scripts"]["spearman_all_vs_all"],
        min_samples = config["correlation"]["min_samples"],
        thresholds  = config["correlation"]["thresholds"],
        p_threshold = config["correlation"]["p_threshold"],
        chunk_size  = config["correlation"]["chunk_size_all_vs_all"],
    log:
        out = "log/05_esvirtu_correlation_analysis/2_all_vs_all/{level}.log",
        err = "log/05_esvirtu_correlation_analysis/2_all_vs_all/{level}.err",
    shell:
        """
        mkdir -p {_OUTDIR_ALL2ALL} log/05_esvirtu_correlation_analysis/2_all_vs_all
        python {params.script} \
            --input {input.table} \
            --output-filtered {output.filt_corr} \
            --min-samples {params.min_samples} \
            --thresholds {params.thresholds} \
            --p-threshold {params.p_threshold} \
            --chunk-size {params.chunk_size} \
            --threads {threads} \
            > {log.out} 2> {log.err}
        """


# ------------------------------------------------------------------ #
# Split filtered Spearman results into known-known / known-unknown   #
# pairs at each r threshold (subspecies and species only)            #
# ------------------------------------------------------------------ #

rule filter_spearman_pairs:
    input:
        filtered = f"{_OUTDIR_ALL2ALL}/{{level}}.all_vs_all.spearman.filtered.tsv",
    output:
        expand(
            f"{_OUTDIR_ANALYSIS}/{{{{level}}}}/{{pair_type}}/r{{threshold}}.tsv",
            pair_type=_PAIR_TYPES,
            threshold=_THRESHOLDS,
        ),
    wildcard_constraints:
        level = "subspecies|species",
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
        script = "scripts/05_esviritu/filter_spearman_pairs.py",
        outdir = f"{_OUTDIR_ANALYSIS}/{{level}}",
    log:
        err = "log/05_esvirtu_correlation_analysis/3_spearman_analysis/{level}.filter.err",
    shell:
        """
        mkdir -p {params.outdir} log/05_esvirtu_correlation_analysis/3_spearman_analysis
        python {params.script} \
            --input  {input.filtered} \
            --outdir {params.outdir} \
            2> {log.err}
        """


# ------------------------------------------------------------------ #
# Annotate known-unknown pairs with RdRP category                    #
# ------------------------------------------------------------------ #

rule annotate_rdrp_category:
    input:
        ku_tsv     = f"{_OUTDIR_ANALYSIS}/{{level}}/known_unknown_pair/r{{threshold}}.tsv",
        tier_tsv   = ancient(_TIER_SUMMARY_TSV),
    output:
        annotated = f"{_OUTDIR_ANNOTATED}/{{level}}/known_unknown_pair/r{{threshold}}.annotated.tsv",
    wildcard_constraints:
        level     = "subspecies|species",
        threshold = "0\\.6|0\\.7|0\\.8|0\\.9",
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = 1,
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    params:
        script = "scripts/05_esviritu/annotate_rdrp_category.py",
    log:
        err = "log/05_esvirtu_correlation_analysis/4_spearman_analysis_RdRp_annotated/{level}.annotate_r{threshold}.err",
    shell:
        """
        mkdir -p $(dirname {output.annotated}) log/05_esvirtu_correlation_analysis/4_spearman_analysis_RdRp_annotated
        python {params.script} \
            --input     {input.ku_tsv} \
            --output    {output.annotated} \
            --tier-tsv  {input.tier_tsv} \
            2> {log.err}
        """


# ------------------------------------------------------------------ #
# FastSpar (SparCC) pipeline — commented out                         #
# ------------------------------------------------------------------ #

# rule fastspar_prepare:
#     input:
#         table = lambda wc: ancient(_RC_INPUT[wc.level]),
#     output:
#         otu = f"{_OUTDIR_SPARCC}/{{level}}/otu_table.tsv",
#     wildcard_constraints:
#         level = "assembly|subspecies|species",
#     conda:
#         "../envs/python.yaml"
#     resources:
#         mem_mb_per_cpu  = config["fastspar"]["memory"],
#         runtime         = 30,
#         slurm_partition = config["fastspar"]["partition"],
#         slurm_account   = config["fastspar"]["account"],
#     params:
#         script             = config["scripts"]["prepare_fastspar_input"],
#         min_samples_others = config["fastspar"]["min_samples_others"],
#     log:
#         err = "log/05_esvirtu_correlation_analysis/3_fastspar/{level}.prepare.err",
#     shell:
#         """
#         mkdir -p {_OUTDIR_SPARCC}/{wildcards.level} log/05_esvirtu_correlation_analysis/3_fastspar
#         python {params.script} \
#             --input {input.table} \
#             --output {output.otu} \
#             --min-samples-others {params.min_samples_others} \
#             2> {log.err}
#         """
#
#
# rule fastspar_run:
#     input:
#         otu = f"{_OUTDIR_SPARCC}/{{level}}/otu_table.tsv",
#     output:
#         corr = f"{_OUTDIR_SPARCC}/{{level}}/correlation.tsv",
#         cov  = f"{_OUTDIR_SPARCC}/{{level}}/covariance.tsv",
#     wildcard_constraints:
#         level = "assembly|subspecies|species",
#     threads: config["fastspar"]["threads"]
#     conda:
#         "../envs/fastspar.yaml"
#     resources:
#         mem_mb_per_cpu  = config["fastspar"]["memory"],
#         runtime         = config["fastspar"]["runtime"],
#         cpus_per_task   = config["fastspar"]["threads"],
#         slurm_partition = config["fastspar"]["partition"],
#         slurm_account   = config["fastspar"]["account"],
#     params:
#         iterations         = config["fastspar"]["iterations"],
#         exclude_iterations = config["fastspar"]["exclude_iterations"],
#         threshold          = config["fastspar"]["threshold"],
#     log:
#         err = "log/05_esvirtu_correlation_analysis/3_fastspar/{level}.fastspar.err",
#     shell:
#         """
#         fastspar \
#             --otu_table {input.otu} \
#             --correlation {output.corr} \
#             --covariance {output.cov} \
#             --iterations {params.iterations} \
#             --exclude_iterations {params.exclude_iterations} \
#             --threshold {params.threshold} \
#             --threads {threads} \
#             2> {log.err}
#         """
#
#
# rule fastspar_bootstrap:
#     input:
#         otu = f"{_OUTDIR_SPARCC}/{{level}}/otu_table.tsv",
#     output:
#         done = f"{_OUTDIR_SPARCC}/{{level}}/bootstrap/.done",
#     wildcard_constraints:
#         level = "assembly|subspecies|species",
#     threads: config["fastspar"]["threads"]
#     conda:
#         "../envs/fastspar.yaml"
#     resources:
#         mem_mb_per_cpu  = config["fastspar"]["memory"],
#         runtime         = config["fastspar"]["runtime"],
#         cpus_per_task   = config["fastspar"]["threads"],
#         slurm_partition = config["fastspar"]["partition"],
#         slurm_account   = config["fastspar"]["account"],
#     params:
#         bootstraps         = config["fastspar"]["bootstraps"],
#         bootstrap_dir      = f"{_OUTDIR_SPARCC}/{{level}}/bootstrap",
#         iterations         = config["fastspar"]["iterations"],
#         exclude_iterations = config["fastspar"]["exclude_iterations"],
#         threshold          = config["fastspar"]["threshold"],
#     log:
#         err = "log/05_esvirtu_correlation_analysis/3_fastspar/{level}.bootstrap.err",
#     shell:
#         """
#         mkdir -p {params.bootstrap_dir}
#         fastspar_bootstrap \
#             --otu_table {input.otu} \
#             --number {params.bootstraps} \
#             --prefix {params.bootstrap_dir}/bootstrap \
#             2> {log.err}
#         parallel -j {threads} \
#             fastspar \
#                 --otu_table {{}} \
#                 --correlation {params.bootstrap_dir}/cor_{{/}} \
#                 --covariance  {params.bootstrap_dir}/cov_{{/}} \
#                 --iterations {params.iterations} \
#                 --exclude_iterations {params.exclude_iterations} \
#                 --threshold {params.threshold} \
#                 --threads 1 \
#             ::: {params.bootstrap_dir}/bootstrap_*.tsv \
#             2>> {log.err}
#         touch {output.done}
#         """
#
#
# rule fastspar_pvalues:
#     input:
#         otu  = f"{_OUTDIR_SPARCC}/{{level}}/otu_table.tsv",
#         corr = f"{_OUTDIR_SPARCC}/{{level}}/correlation.tsv",
#         done = f"{_OUTDIR_SPARCC}/{{level}}/bootstrap/.done",
#     output:
#         pval = f"{_OUTDIR_SPARCC}/{{level}}/pvalues.tsv",
#     wildcard_constraints:
#         level = "assembly|subspecies|species",
#     threads: config["fastspar"]["threads"]
#     conda:
#         "../envs/fastspar.yaml"
#     resources:
#         mem_mb_per_cpu  = config["fastspar"]["memory"],
#         runtime         = config["fastspar"]["runtime"],
#         cpus_per_task   = config["fastspar"]["threads"],
#         slurm_partition = config["fastspar"]["partition"],
#         slurm_account   = config["fastspar"]["account"],
#     params:
#         bootstraps    = config["fastspar"]["bootstraps"],
#         bootstrap_dir = f"{_OUTDIR_SPARCC}/{{level}}/bootstrap",
#     log:
#         err = "log/05_esvirtu_correlation_analysis/3_fastspar/{level}.pvalues.err",
#     shell:
#         """
#         fastspar_pvalues \
#             --otu_table   {input.otu} \
#             --correlation {input.corr} \
#             --prefix      {params.bootstrap_dir}/cor_bootstrap_ \
#             --permutations {params.bootstraps} \
#             --outfile     {output.pval} \
#             --threads     {threads} \
#             2> {log.err}
#         """
#
#
# rule fastspar_parse:
#     input:
#         corr = f"{_OUTDIR_SPARCC}/{{level}}/correlation.tsv",
#         pval = f"{_OUTDIR_SPARCC}/{{level}}/pvalues.tsv",
#     output:
#         filtered = f"{_OUTDIR_SPARCC}/{{level}}/{{level}}.sparcc.filtered.tsv",
#     wildcard_constraints:
#         level = "assembly|subspecies|species",
#     conda:
#         "../envs/python.yaml"
#     resources:
#         mem_mb_per_cpu  = config["fastspar"]["memory"],
#         runtime         = 60,
#         slurm_partition = config["fastspar"]["partition"],
#         slurm_account   = config["fastspar"]["account"],
#     params:
#         script      = config["scripts"]["parse_fastspar_output"],
#         thresholds  = config["fastspar"]["r_threshold"],
#         r_min       = config["fastspar"]["r_min"],
#         p_threshold = config["fastspar"]["p_threshold"],
#     log:
#         err = "log/05_esvirtu_correlation_analysis/3_fastspar/{level}.parse.err",
#     shell:
#         """
#         python {params.script} \
#             --correlation {input.corr} \
#             --pvalues     {input.pval} \
#             --output      {output.filtered} \
#             --thresholds  {params.thresholds} \
#             --r-threshold {params.r_min} \
#             --p-threshold {params.p_threshold} \
#             2> {log.err}
#         """
