configfile: "config/config.yaml"

_FIG_JOB = config["figure_job"]

_RDRP_TAX_PIE_RANKS = list(config["rdrp_taxonomy_pie"]["ranks"])

_TAX_STACK_LEVELS = list(config["rna_virus_taxonomy_stack"]["levels"])
_TAX_STACK_TOP_N  = config["rna_virus_taxonomy_stack"]["top_n"]


rule all:
    input:
        expand(
            "result/100_figure/8_rdrp_taxonomy_pie/rdrp_taxonomy_pie.{rank}.svg",
            rank=_RDRP_TAX_PIE_RANKS,
        ) if config["rdrp_taxonomy_pie"]["run"] else [],
        [
            "result/100_figure/2_contig_tag_barplot/contig_tag_barplot.with_unknown.svg",
            "result/100_figure/2_contig_tag_barplot/contig_tag_barplot.no_unknown.svg",
        ] if config["contig_tag_barplot"]["run"] else [],
        expand(
            "result/100_figure/9_rna_virus_taxonomy_stack/rna_virus_taxonomy_stack.{level}.svg",
            level=_TAX_STACK_LEVELS,
        ) if config["rna_virus_taxonomy_stack"]["run"] else [],


rule plot_rdrp_taxonomy_pie:
    input:
        table = ancient("result/99_tables/3_RDRP_Summary/nr_filtered_protein_contig_cluster.tsv"),
    output:
        svg = "result/100_figure/8_rdrp_taxonomy_pie/rdrp_taxonomy_pie.{rank}.svg",
        tsv = "result/100_figure/8_rdrp_taxonomy_pie/rdrp_taxonomy_pie.{rank}.tsv",
    wildcard_constraints:
        rank = "|".join(_RDRP_TAX_PIE_RANKS),
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
        outdir = "result/100_figure/8_rdrp_taxonomy_pie",
        script = "scripts/100_figures/plot_rdrp_taxonomy_pie.py",
        column = lambda wc: f"rvmt_{wc.rank}",
        title  = lambda wc: f"RDRP DATASET ({wc.rank})",
    log:
        err = "log/100_figure/8_rdrp_taxonomy_pie/{rank}.err",
    shell:
        """
        mkdir -p {params.outdir} log/100_figure/8_rdrp_taxonomy_pie
        python {params.script} \
            --input  {input.table} \
            --column {params.column} \
            --out    {output.svg} \
            --title  "{params.title}" \
            2> {log.err}
        """


rule plot_contig_tag_barplot:
    input:
        contig_table  = ancient("result/99_tables/4_contig_tag_abundance/rna_virus_contig_table.tsv"),
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


rule plot_rna_virus_taxonomy_stack:
    input:
        table = ancient("result/99_tables/5_Recalculated_RPKM_result/1_taxonmy/rna_virus_recalc_tpm.annotated.tsv"),
    output:
        svg   = "result/100_figure/9_rna_virus_taxonomy_stack/rna_virus_taxonomy_stack.{level}.svg",
        table = "result/100_figure/9_rna_virus_taxonomy_stack/rna_virus_taxonomy_stack.{level}.csv",
    wildcard_constraints:
        level = "|".join(_TAX_STACK_LEVELS),
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
        outdir     = "result/100_figure/9_rna_virus_taxonomy_stack",
        script     = "scripts/100_figures/plot_rna_virus_taxonomy_stack.py",
        tax_column = lambda wc: f"tax_{wc.level}",
        top_n      = _TAX_STACK_TOP_N,
    log:
        err = "log/100_figure/9_rna_virus_taxonomy_stack/{level}.err",
    shell:
        """
        mkdir -p {params.outdir} log/100_figure/9_rna_virus_taxonomy_stack
        python {params.script} \
            --input        {input.table} \
            --tax-column   {params.tax_column} \
            --top-n        {params.top_n} \
            --output       {output.svg} \
            --output-table {output.table} \
            2> {log.err}
        """
