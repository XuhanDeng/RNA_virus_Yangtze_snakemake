configfile: "config/config.yaml"

_KU_BASE = "result/05_esvirtu_correlation_analysis/4_spearman_analysis_RdRp_annotated"
_KK_BASE = "result/05_esvirtu_correlation_analysis/3_spearman_analysis"
_OUTDIR  = "result/100_figure/1_spearman_network_figure"

_LEVELS        = config["network_figure"]["levels"]
_THRESHOLDS    = config["network_figure"]["thresholds"]
_VARIANTS      = config["network_figure"]["variants"]
_VARIANT_NAMES = list(_VARIANTS.keys())
_SCRIPT        = "scripts/100_figures/plot_network.py"
_RUN_NETWORK   = config["network_figure"]["run"]

_CONTIG_TABLE       = "result/99_tables/3_RNA_virus_table/rna_virus_contig_table.tsv"
_ESVIRTU_TABLE      = "result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv"
_TAG_BARPLOT_OUTDIR = "result/100_figure/2_contig_tag_barplot"
_TAG_BARPLOT_SCRIPT = "scripts/100_figures/plot_contig_tag_barplot.py"
_RUN_TAG_BARPLOT    = config.get("contig_tag_barplot", {}).get("run", False)

# ── host lineage heatmap (5_Recalculated_RPKM_result: full / restricted / tier1only) ──
_RECALC_RPKM_DIR = "result/99_tables/5_Recalculated_RPKM_result"
_HOST_LINEAGE_OUTDIR = "result/100_figure/3_host_lineage_heatmap"
_HOST_LINEAGE_SCRIPT = "scripts/100_figures/plot_host_lineage_barplot.py"

_HOST_LINEAGE_CFG    = config.get("host_lineage_barplot", {})
_RUN_HOST_LINEAGE    = _HOST_LINEAGE_CFG.get("run", False)
_HOST_LINEAGE_TOP_N  = _HOST_LINEAGE_CFG.get("top_n", 15)
_HOST_LINEAGE_VALUE_TYPE = _HOST_LINEAGE_CFG.get("value_type", "both")  # rpkmf | tpm | both

_HOST_LINEAGE_VALUE_TYPES = (
    ["rpkmf", "tpm"] if _HOST_LINEAGE_VALUE_TYPE == "both" else [_HOST_LINEAGE_VALUE_TYPE]
)

# variant name -> subdirectory under 5_Recalculated_RPKM_result
_HOST_LINEAGE_VARIANT_DIRS = {"full": "1_RNA_virus_Recalculated_RPKM"}
if config.get("recalc_rpkmf", {}).get("run_restricted", False):
    _HOST_LINEAGE_VARIANT_DIRS["restricted"] = "2_RNA_virus_Recalculated_RPKM_restricted"
if config.get("recalc_rpkmf", {}).get("run_tier1_only", False):
    _HOST_LINEAGE_VARIANT_DIRS["tier1only"] = "3_RNA_virus_Recalculated_RPKM_tier1only"
_HOST_LINEAGE_VARIANTS = list(_HOST_LINEAGE_VARIANT_DIRS.keys())

# ── known-virus heatmap (5_Recalculated_RPKM_result: full / restricted / tier1only) ──
_KNOWN_VIRUS_OUTDIR = "result/100_figure/4_known_virus_heatmap"
_KNOWN_VIRUS_SCRIPT  = "scripts/100_figures/plot_known_virus_heatmap.py"

_KNOWN_VIRUS_CFG   = config.get("known_virus_heatmap", {})
_RUN_KNOWN_VIRUS   = _KNOWN_VIRUS_CFG.get("run", False)
_KNOWN_VIRUS_TOP_N = _KNOWN_VIRUS_CFG.get("top_n", 15)
_KNOWN_VIRUS_VALUE_TYPE = _KNOWN_VIRUS_CFG.get("value_type", "both")  # rpkmf | tpm | both

_KNOWN_VIRUS_VALUE_TYPES = (
    ["rpkmf", "tpm"] if _KNOWN_VIRUS_VALUE_TYPE == "both" else [_KNOWN_VIRUS_VALUE_TYPE]
)

# reuse the same variant -> subdirectory mapping as the host lineage heatmap
_KNOWN_VIRUS_VARIANT_DIRS = _HOST_LINEAGE_VARIANT_DIRS
_KNOWN_VIRUS_VARIANTS = list(_KNOWN_VIRUS_VARIANT_DIRS.keys())


rule all:
    input:
        expand(f"{_OUTDIR}/{{level}}/{{variant}}/ku{{threshold}}_with_label.svg",
               level=_LEVELS, variant=_VARIANT_NAMES, threshold=_THRESHOLDS) if _RUN_NETWORK else [],
        expand(f"{_OUTDIR}/{{level}}/{{variant}}/ku{{threshold}}_no_label.svg",
               level=_LEVELS, variant=_VARIANT_NAMES, threshold=_THRESHOLDS) if _RUN_NETWORK else [],
        expand(f"{_OUTDIR}/{{level}}/{{variant}}/ku{{threshold}}_with_label.html",
               level=_LEVELS, variant=_VARIANT_NAMES, threshold=_THRESHOLDS) if _RUN_NETWORK else [],
        expand(f"{_OUTDIR}/{{level}}/{{variant}}/ku{{threshold}}_no_label.html",
               level=_LEVELS, variant=_VARIANT_NAMES, threshold=_THRESHOLDS) if _RUN_NETWORK else [],
        [
            f"{_TAG_BARPLOT_OUTDIR}/contig_tag_barplot.with_unknown.svg",
            f"{_TAG_BARPLOT_OUTDIR}/contig_tag_barplot.no_unknown.svg",
        ] if _RUN_TAG_BARPLOT else [],
        expand(
            f"{_HOST_LINEAGE_OUTDIR}/{{variant}}/host_lineage_heatmap.{{value_type}}.svg",
            variant=_HOST_LINEAGE_VARIANTS,
            value_type=_HOST_LINEAGE_VALUE_TYPES,
        ) if _RUN_HOST_LINEAGE else [],
        expand(
            f"{_KNOWN_VIRUS_OUTDIR}/{{variant}}/known_virus_heatmap.{{value_type}}.svg",
            variant=_KNOWN_VIRUS_VARIANTS,
            value_type=_KNOWN_VIRUS_VALUE_TYPES,
        ) if _RUN_KNOWN_VIRUS else [],


rule plot_network:
    input:
        ku_dir = ancient(f"{_KU_BASE}/{{level}}/known_unknown_pair"),
        kk_dir = ancient(f"{_KK_BASE}/{{level}}/known_known_pair"),
    output:
        svg_label    = f"{_OUTDIR}/{{level}}/{{variant}}/ku{{threshold}}_with_label.svg",
        svg_nolabel  = f"{_OUTDIR}/{{level}}/{{variant}}/ku{{threshold}}_no_label.svg",
        html_label   = f"{_OUTDIR}/{{level}}/{{variant}}/ku{{threshold}}_with_label.html",
        html_nolabel = f"{_OUTDIR}/{{level}}/{{variant}}/ku{{threshold}}_no_label.html",
    wildcard_constraints:
        level     = "subspecies|species",
        threshold = "0\\.6|0\\.7|0\\.8|0\\.9",
        variant   = "|".join(_VARIANT_NAMES),
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
        kk_threshold     = lambda wc: _VARIANTS[wc.variant]["kk_threshold"],
        hide_not_in_list = lambda wc: _VARIANTS[wc.variant]["hide_not_in_list"],
        kk_only          = lambda wc: _VARIANTS[wc.variant]["kk_only"],
    log:
        err = "log/100_network_figure/1_spearman_network_figure/{level}/{variant}/ku{threshold}.err",
    shell:
        """
        mkdir -p {_OUTDIR}/{wildcards.level}/{wildcards.variant} \
                 log/100_network_figure/1_spearman_network_figure/{wildcards.level}/{wildcards.variant}

        kk_arg=""
        if [ "{params.kk_threshold}" = "all" ]; then
            kk_arg="--all-kk"
        elif [ -n "{params.kk_threshold}" ]; then
            kk_arg="--kk-threshold {params.kk_threshold}"
        fi

        nil_arg=""
        if [ "{params.hide_not_in_list}" = "True" ]; then
            nil_arg="--hide-not-in-list"
        fi

        kk_only_arg=""
        if [ "{params.kk_only}" = "True" ]; then
            kk_only_arg="--kk-only"
        fi

        python {_SCRIPT} \
            --ku-dir        {input.ku_dir} \
            --kk-dir        {input.kk_dir} \
            --min-threshold {wildcards.threshold} \
            $kk_arg \
            $nil_arg \
            $kk_only_arg \
            --out-svg-label    {output.svg_label} \
            --out-svg-nolabel  {output.svg_nolabel} \
            --out-html-label   {output.html_label} \
            --out-html-nolabel {output.html_nolabel} \
            2> {log.err}
        """


rule plot_contig_tag_barplot:
    input:
        contig_table  = ancient(_CONTIG_TABLE),
        esvirtu_table = ancient(_ESVIRTU_TABLE),
    output:
        with_unknown       = f"{_TAG_BARPLOT_OUTDIR}/contig_tag_barplot.with_unknown.svg",
        no_unknown         = f"{_TAG_BARPLOT_OUTDIR}/contig_tag_barplot.no_unknown.svg",
        table_with_unknown = f"{_TAG_BARPLOT_OUTDIR}/contig_tag_barplot.with_unknown.csv",
        table_no_unknown   = f"{_TAG_BARPLOT_OUTDIR}/contig_tag_barplot.no_unknown.csv",
    conda:
        "../envs/python.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        err = "log/100_figure/2_contig_tag_barplot.err",
    shell:
        """
        mkdir -p {_TAG_BARPLOT_OUTDIR} log/100_figure

        python {_TAG_BARPLOT_SCRIPT} \
            --contig-table  {input.contig_table} \
            --esvirtu-table {input.esvirtu_table} \
            --output        {output.with_unknown} \
            --output-table  {output.table_with_unknown} \
            --include-unknown \
            2> {log.err}

        python {_TAG_BARPLOT_SCRIPT} \
            --contig-table  {input.contig_table} \
            --esvirtu-table {input.esvirtu_table} \
            --output        {output.no_unknown} \
            --output-table  {output.table_no_unknown} \
            2>> {log.err}
        """


rule plot_host_lineage_heatmap:
    input:
        table = lambda wc: ancient(
            f"{_RECALC_RPKM_DIR}/{_HOST_LINEAGE_VARIANT_DIRS[wc.variant]}/rna_virus_recalc_{wc.value_type}.tsv"
        ),
    output:
        svg   = f"{_HOST_LINEAGE_OUTDIR}/{{variant}}/host_lineage_heatmap.{{value_type}}.svg",
        table = f"{_HOST_LINEAGE_OUTDIR}/{{variant}}/host_lineage_heatmap.{{value_type}}.all_hosts.csv",
    wildcard_constraints:
        variant    = "|".join(_HOST_LINEAGE_VARIANTS) if _HOST_LINEAGE_VARIANTS else "full",
        value_type = "rpkmf|tpm",
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
        value_suffix = lambda wc: f"_{wc.value_type}",
        top_n        = _HOST_LINEAGE_TOP_N,
    log:
        err = "log/100_figure/3_host_lineage_heatmap/{variant}.{value_type}.err",
    shell:
        """
        mkdir -p {_HOST_LINEAGE_OUTDIR}/{wildcards.variant} log/100_figure/3_host_lineage_heatmap

        python {_HOST_LINEAGE_SCRIPT} \
            --input         {input.table} \
            --value-suffix  {params.value_suffix} \
            --top-n         {params.top_n} \
            --output        {output.svg} \
            --output-table  {output.table} \
            2> {log.err}
        """


rule plot_known_virus_heatmap:
    input:
        table = lambda wc: ancient(
            f"{_RECALC_RPKM_DIR}/{_KNOWN_VIRUS_VARIANT_DIRS[wc.variant]}/rna_virus_recalc_{wc.value_type}.tsv"
        ),
    output:
        svg = f"{_KNOWN_VIRUS_OUTDIR}/{{variant}}/known_virus_heatmap.{{value_type}}.svg",
    wildcard_constraints:
        variant    = "|".join(_KNOWN_VIRUS_VARIANTS) if _KNOWN_VIRUS_VARIANTS else "full",
        value_type = "rpkmf|tpm",
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
        value_suffix = lambda wc: f"_{wc.value_type}",
        top_n        = _KNOWN_VIRUS_TOP_N,
    log:
        err = "log/100_figure/4_known_virus_heatmap/{variant}.{value_type}.err",
    shell:
        """
        mkdir -p {_KNOWN_VIRUS_OUTDIR}/{wildcards.variant} log/100_figure/4_known_virus_heatmap

        python {_KNOWN_VIRUS_SCRIPT} \
            --input         {input.table} \
            --value-suffix  {params.value_suffix} \
            --top-n         {params.top_n} \
            --output        {output.svg} \
            2> {log.err}
        """
