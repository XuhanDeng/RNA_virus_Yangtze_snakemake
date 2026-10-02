configfile: "config/config.yaml"
import re

# Only kept here because expand() / wildcard_constraints need them at parse time
_NETWORK_THRESHOLDS  = config["network_figure"]["thresholds"]
_NETWORK_VARIANTS    = list(config["network_figure"]["variants"].keys())
_NETWORK_VALUE_TYPES = config["network_figure"].get("value_types", ["rpkmf"])

_HOST_LINEAGE_VALUE_TYPE  = config["host_lineage_barplot"]["value_type"]
_HOST_LINEAGE_VALUE_TYPES = ["rpkmf", "tpm"] if _HOST_LINEAGE_VALUE_TYPE == "both" else [_HOST_LINEAGE_VALUE_TYPE]

_KNOWN_VIRUS_VALUE_TYPE  = config["known_virus_heatmap"]["value_type"]
_KNOWN_VIRUS_VALUE_TYPES = ["rpkmf", "tpm"] if _KNOWN_VIRUS_VALUE_TYPE == "both" else [_KNOWN_VIRUS_VALUE_TYPE]

_COABUNDANCE_THRESHOLDS = [str(t) for t in config["coabundance_figure"]["thresholds"]]

_FIG_JOB = config["figure_job"]

_CORR_BASE     = "result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation"
_CORR_BASE_TPM = "result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation_TPM"
_CORR_KU = {
    "rpkmf": _CORR_BASE     + "/3_spearman_analysis_RdRp_annotated",
    "tpm":   _CORR_BASE_TPM + "/3_spearman_analysis_RdRp_annotated",
}
_CORR_KK = {
    "rpkmf": _CORR_BASE     + "/2_spearman_analysis",
    "tpm":   _CORR_BASE_TPM + "/2_spearman_analysis",
}


rule all:
    input:
        expand(
            "result/100_figure/1_spearman_network_figure/{value_type}/{variant}/ku{threshold}_with_label.svg",
            value_type=_NETWORK_VALUE_TYPES, variant=_NETWORK_VARIANTS, threshold=_NETWORK_THRESHOLDS,
        ) if config["network_figure"]["run"] else [],
        expand(
            "result/100_figure/1_spearman_network_figure/{value_type}/{variant}/ku{threshold}_no_label.svg",
            value_type=_NETWORK_VALUE_TYPES, variant=_NETWORK_VARIANTS, threshold=_NETWORK_THRESHOLDS,
        ) if config["network_figure"]["run"] else [],
        expand(
            "result/100_figure/1_spearman_network_figure/{value_type}/{variant}/ku{threshold}_with_label.html",
            value_type=_NETWORK_VALUE_TYPES, variant=_NETWORK_VARIANTS, threshold=_NETWORK_THRESHOLDS,
        ) if config["network_figure"]["run"] else [],
        expand(
            "result/100_figure/1_spearman_network_figure/{value_type}/{variant}/ku{threshold}_no_label.html",
            value_type=_NETWORK_VALUE_TYPES, variant=_NETWORK_VARIANTS, threshold=_NETWORK_THRESHOLDS,
        ) if config["network_figure"]["run"] else [],
        [
            "result/100_figure/2_contig_tag_barplot/contig_tag_barplot.with_unknown.svg",
            "result/100_figure/2_contig_tag_barplot/contig_tag_barplot.no_unknown.svg",
        ] if config["contig_tag_barplot"]["run"] else [],
        expand(
            "result/100_figure/3_host_lineage_heatmap/host_lineage_heatmap.{value_type}.svg",
            value_type=_HOST_LINEAGE_VALUE_TYPES,
        ) if config["host_lineage_barplot"]["run"] else [],
        expand(
            "result/100_figure/4_known_virus_heatmap/known_virus_heatmap.{value_type}.svg",
            value_type=_KNOWN_VIRUS_VALUE_TYPES,
        ) if config["known_virus_heatmap"]["run"] else [],
        [
            "result/100_figure/5_contamination/contamination.svg",
            "result/100_figure/5_contamination/contamination.tsv",
        ] if config["contamination_figure"]["run"] else [],
        expand(
            "result/100_figure/6_gappa_phylum/{batch}/phylum_barplot.svg",
            batch=config["gappa_phylum_figure"]["batches"],
        ) if config["gappa_phylum_figure"]["run"] else [],
        expand(
            "result/100_figure/7_coabundance/{threshold}/.done",
            threshold=_COABUNDANCE_THRESHOLDS,
        ) if config["coabundance_figure"]["run"] else [],


rule plot_network:
    input:
        ku_dir = lambda wc: ancient(f"{_CORR_KU[wc.value_type]}/known_unknown_pair"),
        kk_dir = lambda wc: ancient(f"{_CORR_KK[wc.value_type]}/known_known_pair"),
    output:
        svg_label    = "result/100_figure/1_spearman_network_figure/{value_type}/{variant}/ku{threshold}_with_label.svg",
        svg_nolabel  = "result/100_figure/1_spearman_network_figure/{value_type}/{variant}/ku{threshold}_no_label.svg",
        html_label   = "result/100_figure/1_spearman_network_figure/{value_type}/{variant}/ku{threshold}_with_label.html",
        html_nolabel = "result/100_figure/1_spearman_network_figure/{value_type}/{variant}/ku{threshold}_no_label.html",
    wildcard_constraints:
        value_type = "rpkmf|tpm",
        threshold  = "0\\.6|0\\.7|0\\.8|0\\.9",
        variant    = "|".join(_NETWORK_VARIANTS),
    conda:
        "../envs/python.yaml"
    threads: _FIG_JOB["threads"]
    resources:
        mem_mb_per_cpu  = _FIG_JOB["memory"],
        runtime         = _FIG_JOB["runtime"],
        cpus_per_task   = _FIG_JOB["threads"],
        slurm_partition = _FIG_JOB["partition"],
        slurm_account   = _FIG_JOB["account"],
    params:
        outdir           = "result/100_figure/1_spearman_network_figure",
        script           = "scripts/100_figures/plot_network.py",
        kk_threshold     = lambda wc: config["network_figure"]["variants"][wc.variant]["kk_threshold"],
        hide_not_in_list = lambda wc: config["network_figure"]["variants"][wc.variant]["hide_not_in_list"],
        kk_only          = lambda wc: config["network_figure"]["variants"][wc.variant]["kk_only"],
    log:
        err = "log/100_figure/1_spearman_network_figure/{value_type}/{variant}/ku{threshold}.err",
    shell:
        """
        mkdir -p {params.outdir}/{wildcards.value_type}/{wildcards.variant} \
                 log/100_figure/1_spearman_network_figure/{wildcards.value_type}/{wildcards.variant}

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

        python {params.script} \
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
        contig_table  = ancient("result/99_tables/3_RNA_virus_table/rna_virus_contig_table.tsv"),
        esvirtu_table = ancient("result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv"),
    output:
        with_unknown       = "result/100_figure/2_contig_tag_barplot/contig_tag_barplot.with_unknown.svg",
        no_unknown         = "result/100_figure/2_contig_tag_barplot/contig_tag_barplot.no_unknown.svg",
        table_with_unknown = "result/100_figure/2_contig_tag_barplot/contig_tag_barplot.with_unknown.csv",
        table_no_unknown   = "result/100_figure/2_contig_tag_barplot/contig_tag_barplot.no_unknown.csv",
    conda:
        "../envs/python.yaml"
    threads: _FIG_JOB["threads"]
    resources:
        mem_mb_per_cpu  = _FIG_JOB["memory"],
        runtime         = _FIG_JOB["runtime"],
        cpus_per_task   = _FIG_JOB["threads"],
        slurm_partition = _FIG_JOB["partition"],
        slurm_account   = _FIG_JOB["account"],
    params:
        outdir = "result/100_figure/2_contig_tag_barplot",
        script = "scripts/100_figures/plot_contig_tag_barplot.py",
    log:
        err = "log/100_figure/2_contig_tag_barplot.err",
    shell:
        """
        mkdir -p {params.outdir} log/100_figure

        python {params.script} \
            --contig-table  {input.contig_table} \
            --esvirtu-table {input.esvirtu_table} \
            --output        {output.with_unknown} \
            --output-table  {output.table_with_unknown} \
            --include-unknown \
            2> {log.err}

        python {params.script} \
            --contig-table  {input.contig_table} \
            --esvirtu-table {input.esvirtu_table} \
            --output        {output.no_unknown} \
            --output-table  {output.table_no_unknown} \
            2>> {log.err}
        """


rule plot_host_lineage_heatmap:
    input:
        table = lambda wc: ancient(
            f"result/99_tables/5_Recalculated_RPKM_result/rna_virus_recalc_{wc.value_type}.tsv"
        ),
    output:
        svg   = "result/100_figure/3_host_lineage_heatmap/host_lineage_heatmap.{value_type}.svg",
        table = "result/100_figure/3_host_lineage_heatmap/host_lineage_heatmap.{value_type}.all_hosts.csv",
    wildcard_constraints:
        value_type = "rpkmf|tpm",
    conda:
        "../envs/python.yaml"
    threads: _FIG_JOB["threads"]
    resources:
        mem_mb_per_cpu  = _FIG_JOB["memory"],
        runtime         = _FIG_JOB["runtime"],
        cpus_per_task   = _FIG_JOB["threads"],
        slurm_partition = _FIG_JOB["partition"],
        slurm_account   = _FIG_JOB["account"],
    params:
        outdir       = "result/100_figure/3_host_lineage_heatmap",
        script       = "scripts/100_figures/plot_host_lineage_barplot.py",
        value_suffix = lambda wc: f"_{wc.value_type}",
        top_n        = config["host_lineage_barplot"]["top_n"],
    log:
        err = "log/100_figure/3_host_lineage_heatmap/{value_type}.err",
    shell:
        """
        mkdir -p {params.outdir} log/100_figure/3_host_lineage_heatmap

        python {params.script} \
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
            f"result/99_tables/5_Recalculated_RPKM_result/rna_virus_recalc_{wc.value_type}.tsv"
        ),
    output:
        svg = "result/100_figure/4_known_virus_heatmap/known_virus_heatmap.{value_type}.svg",
    wildcard_constraints:
        value_type = "rpkmf|tpm",
    conda:
        "../envs/python.yaml"
    threads: _FIG_JOB["threads"]
    resources:
        mem_mb_per_cpu  = _FIG_JOB["memory"],
        runtime         = _FIG_JOB["runtime"],
        cpus_per_task   = _FIG_JOB["threads"],
        slurm_partition = _FIG_JOB["partition"],
        slurm_account   = _FIG_JOB["account"],
    params:
        outdir       = "result/100_figure/4_known_virus_heatmap",
        script       = "scripts/100_figures/plot_known_virus_heatmap.py",
        value_suffix = lambda wc: f"_{wc.value_type}",
        top_n        = config["known_virus_heatmap"]["top_n"],
    log:
        err = "log/100_figure/4_known_virus_heatmap/{value_type}.err",
    shell:
        """
        mkdir -p {params.outdir} log/100_figure/4_known_virus_heatmap

        python {params.script} \
            --input         {input.table} \
            --value-suffix  {params.value_suffix} \
            --top-n         {params.top_n} \
            --output        {output.svg} \
            2> {log.err}
        """


rule plot_contamination:
    input:
        kraken_dir = ancient(config["contamination_figure"]["kraken_dir"]),
    output:
        svg = "result/100_figure/5_contamination/contamination.svg",
        tsv = "result/100_figure/5_contamination/contamination.tsv",
    conda:
        "../envs/python.yaml"
    threads: _FIG_JOB["threads"]
    resources:
        mem_mb_per_cpu  = _FIG_JOB["memory"],
        runtime         = _FIG_JOB["runtime"],
        cpus_per_task   = _FIG_JOB["threads"],
        slurm_partition = _FIG_JOB["partition"],
        slurm_account   = _FIG_JOB["account"],
    params:
        outdir    = "result/100_figure/5_contamination",
        script    = "scripts/100_figures/plot_contamination.py",
        threshold = config["contamination_figure"]["threshold"],
    log:
        err = "log/100_figure/5_contamination/contamination.err",
    shell:
        """
        mkdir -p {params.outdir} log/100_figure/5_contamination
        python {params.script} \
            --kraken-dir {input.kraken_dir} \
            --out-svg    {output.svg} \
            --threshold  {params.threshold} \
            2> {log.err}
        """


rule plot_gappa_phylum:
    input:
        tsv = config["gappa_phylum_figure"]["indir"] + "/{batch}/phylum_annotation.tsv",
    output:
        svg = "result/100_figure/6_gappa_phylum/{batch}/phylum_barplot.svg",
        tsv = "result/100_figure/6_gappa_phylum/{batch}/phylum_barplot.tsv",
    conda:
        "../envs/python.yaml"
    threads: _FIG_JOB["threads"]
    resources:
        mem_mb_per_cpu  = _FIG_JOB["memory"],
        runtime         = _FIG_JOB["runtime"],
        cpus_per_task   = _FIG_JOB["threads"],
        slurm_partition = _FIG_JOB["partition"],
        slurm_account   = _FIG_JOB["account"],
    params:
        outdir = "result/100_figure/6_gappa_phylum",
        script = "scripts/100_figures/plot_gappa_phylum.py",
        title  = lambda wc: f"Phylum Assignment — {wc.batch} (gappa EPA-ng)",
    log:
        err = "log/100_figure/6_gappa_phylum/{batch}.err",
    shell:
        """
        mkdir -p $(dirname {output.svg}) log/100_figure/6_gappa_phylum
        python {params.script} \
            --input {input.tsv} \
            --out   {output.svg} \
            --title "{params.title}" \
            2> {log.err}
        """


rule plot_coabundance:
    input:
        tpm    = ancient(config["coabundance_figure"]["tpm"]),
        kk_dir = ancient(config["coabundance_figure"]["kk_dir"]),
        ku_dir = ancient(config["coabundance_figure"]["ku_dir"]),
    output:
        done = "result/100_figure/7_coabundance/{threshold}/.done",
    wildcard_constraints:
        threshold = "|".join(re.escape(t) for t in _COABUNDANCE_THRESHOLDS),
    conda:
        "../envs/python.yaml"
    threads: _FIG_JOB["threads"]
    resources:
        mem_mb_per_cpu  = _FIG_JOB["memory"],
        runtime         = _FIG_JOB["runtime"],
        cpus_per_task   = _FIG_JOB["threads"],
        slurm_partition = _FIG_JOB["partition"],
        slurm_account   = _FIG_JOB["account"],
    params:
        outdir = lambda wc: f"result/100_figure/7_coabundance/{wc.threshold}",
        script = "scripts/100_figures/plot_virus_coabundance.py",
    log:
        err = "log/100_figure/7_coabundance/{threshold}.err",
    shell:
        """
        mkdir -p {params.outdir} log/100_figure/7_coabundance
        python {params.script} \
            --tpm       {input.tpm} \
            --kk-dir    {input.kk_dir} \
            --ku-dir    {input.ku_dir} \
            --threshold {wildcards.threshold} \
            --outdir    {params.outdir} \
            2> {log.err}
        touch {output.done}
        """
