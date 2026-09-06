configfile: "config/config.yaml"

_SCRATCH = "/scratch/xddeng/yangtze/RNA/my_rna/result/03_RDRP_identification/5_RdRp_contig"
_OUTDIR  = "result/test_lucaprot_palmscan"

_TARGETS = {
    "lucaprot":         f"{_SCRATCH}/lucaprot.fasta",
    "RDRPCatch_lucaprot": f"{_SCRATCH}/RDRPCatch_lucaprot.fasta",
}

rule all:
    input:
        expand(f"{_OUTDIR}/{{name}}_palmscan_hits.tsv", name=_TARGETS.keys()),


rule palm_annot_test:
    input:
        aa = lambda wc: ancient(_TARGETS[wc.name]),
    output:
        tsv     = f"{_OUTDIR}/{{name}}_palmscan_hits.tsv",
        fev     = f"{_OUTDIR}/{{name}}_palm_annot.fev",
        rdrp_aa = f"{_OUTDIR}/{{name}}_rdrp_trimmed.faa",
    wildcard_constraints:
        name = "|".join(_TARGETS.keys()),
    conda:
        "../envs/palm_annot.yaml"
    threads: config["palm_annot"]["threads"]
    resources:
        mem_mb_per_cpu  = config["palm_annot"]["memory"],
        runtime         = config["palm_annot"]["runtime"],
        cpus_per_task   = config["palm_annot"]["threads"],
        slurm_partition = config["palm_annot"]["partition"],
        slurm_account   = config["palm_annot"]["account"],
    params:
        palm_annot_py = config["palm_annot"]["install_dir"] + "/py/palm_annot.py",
        fev2tsv_py    = config["palm_annot"]["install_dir"] + "/py/fev2tsv.py",
        seqtype       = "nt",
        minscore      = config["palm_annot"]["minscore"],
        minpssmscore  = config["palm_annot"]["minpssmscore"],
        tmpdir        = "tmp/palm_annot/test/{name}",
    log:
        out = "log/test_lucaprot_palmscan/{name}.log",
        err = "log/test_lucaprot_palmscan/{name}.err",
    shell:
        """
        mkdir -p {_OUTDIR} {params.tmpdir} log/test_lucaprot_palmscan
        python {params.palm_annot_py} \
            --input        {input.aa} \
            --seqtype      {params.seqtype} \
            --fev          {output.fev} \
            --rdrp         {output.rdrp_aa} \
            --minscore     {params.minscore} \
            --minpssmscore {params.minpssmscore} \
            --threads      {threads} \
            --tmpdir       {params.tmpdir} \
            > {log.out} 2> {log.err}
        python {params.fev2tsv_py} \
            --input  {output.fev} \
            --output {output.tsv} \
            --header yes \
            >> {log.out} 2>> {log.err}
        """
