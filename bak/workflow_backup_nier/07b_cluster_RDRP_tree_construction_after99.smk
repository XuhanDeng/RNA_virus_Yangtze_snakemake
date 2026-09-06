# RdRP phylogenetic PLACEMENT pipeline — clustered contigs only (after workflow 99)
#
# Standalone pipeline; does NOT depend on 07_RDRP_tree_construction.smk.
# Inputs come from workflow 03 (tier FAA files) and workflow 99 (clustered table).
#
# Differences from 07:
#   - Only contigs present in rna_virus_contig_table.tsv (clustered set) are used
#     as queries. The full tier1a/1b/2 FAA files are pre-filtered in Step 0b.
#   - ICTV / ESvirtu reference candidates are still included when enabled in
#     config, because they are not assembly contigs and are not filtered.
#   - Outputs go to result/07b_cluster_RDRP_tree_construction/
#   - Log files go to log/07b_cluster_RDRP_tree_construction/
#
# Run AFTER 99_tables.smk has produced rna_virus_contig_table.tsv.

import os

configfile: "config/config.yaml"

_MOTIF_RESULTS      = "result/03_RDRP_identification/4_motif_search/motif_results"
_TIER1A_FAA         = _MOTIF_RESULTS + "/tier1a_ABCD.faa"
_TIER1B_FAA         = _MOTIF_RESULTS + "/tier1b_CABD_depermuted.faa"
_TIER2_FAA          = _MOTIF_RESULTS + "/tier2_3motif.faa"

_ICTV_MOTIF_RESULTS = "result/03_RDRP_identification/4_motif_search/ICTV/motif_results"
_ICTV_TIER1A_FAA    = _ICTV_MOTIF_RESULTS + "/tier1a_ABCD.faa"
_ICTV_TIER1B_FAA    = _ICTV_MOTIF_RESULTS + "/tier1b_CABD_depermuted.faa"
_ICTV_TIER2_FAA     = _ICTV_MOTIF_RESULTS + "/tier2_3motif.faa"

_ESV_MOTIF_RESULTS  = "result/04_03_esvirtu_rdrp_identification/7_motif_search/motif_results"
_ESV_TIER1A_FAA     = _ESV_MOTIF_RESULTS + "/tier1a_ABCD.faa"
_ESV_TIER1B_FAA     = _ESV_MOTIF_RESULTS + "/tier1b_CABD_depermuted.faa"
_ESV_TIER2_FAA      = _ESV_MOTIF_RESULTS + "/tier2_3motif.faa"

# Clustered contig table produced by workflow 99
_RNA_VIRUS_TABLE    = "result/99_tables/3_RNA_virus_table/rna_virus_contig_table.tsv"

_OUTDIR      = "result/07b_cluster_RDRP_tree_construction"
_REF_MSA     = config["epa_ng"]["ref_msa"]
_REF_TREE    = config["epa_ng"]["ref_tree"]
_MODEL       = config["epa_ng"]["model"]
_RVMT_INFO_TSV = config["epa_ng"]["rvmt_info_tsv"]
_TAXON_FILE  = _OUTDIR + "/0_taxon/RCR90_taxonomy.tsv"

_INCLUDE_ICTV     = config["epa_ng"].get("include_ictv", False)
_INCLUDE_ESVIRTU  = config["epa_ng"].get("include_esvirtu", False)
_INCLUDE_FASTTREE = config["epa_ng"].get("include_fasttree", False)

# Contig FAA files are filtered versions (0b_filtered); ICTV/ESvirtu are staged
# directly from their source directories (not filtered — they are references).
_STAGED_INPUTS = {
    "contig_tier1a_ABCD.faa":            _OUTDIR + "/0b_filtered/contig_tier1a_ABCD.faa",
    "contig_tier1b_CABD_depermuted.faa": _OUTDIR + "/0b_filtered/contig_tier1b_CABD_depermuted.faa",
    "contig_tier2_3motif.faa":           _OUTDIR + "/0b_filtered/contig_tier2_3motif.faa",
}
if _INCLUDE_ICTV:
    _STAGED_INPUTS.update({
        "ictv_tier1a_ABCD.faa":            _ICTV_TIER1A_FAA,
        "ictv_tier1b_CABD_depermuted.faa": _ICTV_TIER1B_FAA,
        "ictv_tier2_3motif.faa":           _ICTV_TIER2_FAA,
    })
if _INCLUDE_ESVIRTU:
    _STAGED_INPUTS.update({
        "esvirtu_tier1a_ABCD.faa":            _ESV_TIER1A_FAA,
        "esvirtu_tier1b_CABD_depermuted.faa": _ESV_TIER1B_FAA,
        "esvirtu_tier2_3motif.faa":           _ESV_TIER2_FAA,
    })

_BATCHES = {
    "tier1ab":  [_OUTDIR + "/0_inputdata/contig_tier1a_ABCD.faa",
                _OUTDIR + "/0_inputdata/contig_tier1b_CABD_depermuted.faa"],
    "tier1ab2": [_OUTDIR + "/0_inputdata/contig_tier1a_ABCD.faa",
                _OUTDIR + "/0_inputdata/contig_tier1b_CABD_depermuted.faa",
                _OUTDIR + "/0_inputdata/contig_tier2_3motif.faa"],
}
if _INCLUDE_ICTV:
    _BATCHES["tier1ab"]  += [_OUTDIR + "/0_inputdata/ictv_tier1a_ABCD.faa",
                             _OUTDIR + "/0_inputdata/ictv_tier1b_CABD_depermuted.faa"]
    _BATCHES["tier1ab2"] += [_OUTDIR + "/0_inputdata/ictv_tier1a_ABCD.faa",
                             _OUTDIR + "/0_inputdata/ictv_tier1b_CABD_depermuted.faa",
                             _OUTDIR + "/0_inputdata/ictv_tier2_3motif.faa"]
if _INCLUDE_ESVIRTU:
    _BATCHES["tier1ab"]  += [_OUTDIR + "/0_inputdata/esvirtu_tier1a_ABCD.faa",
                             _OUTDIR + "/0_inputdata/esvirtu_tier1b_CABD_depermuted.faa"]
    _BATCHES["tier1ab2"] += [_OUTDIR + "/0_inputdata/esvirtu_tier1a_ABCD.faa",
                             _OUTDIR + "/0_inputdata/esvirtu_tier1b_CABD_depermuted.faa",
                             _OUTDIR + "/0_inputdata/esvirtu_tier2_3motif.faa"]


rule all:
    input:
        _TAXON_FILE,
        # Step 0a — RdRp ID list derived from clustered table + tier_summary
        _OUTDIR + "/0b_filtered/clustered_rdrp_ids.txt",
        # Step 0b — filtered contig FAA files
        expand(_OUTDIR + "/0b_filtered/contig_{tier}.faa",
               tier=["tier1a_ABCD", "tier1b_CABD_depermuted", "tier2_3motif"]),
        # Step 0 — staged input copies (filtered contigs + optional ICTV/ESvirtu)
        expand(_OUTDIR + "/0_inputdata/{f}", f=list(_STAGED_INPUTS.keys())),
        # Step 1 — reference HMM
        _OUTDIR + "/1_ref_hmm/RCR90.hmm",
        # Step 2-3 — per-batch aligned + split
        expand(_OUTDIR + "/3_split/{batch}/query.fasta",     batch=list(_BATCHES.keys())),
        expand(_OUTDIR + "/3_split/{batch}/reference.fasta", batch=list(_BATCHES.keys())),
        # Step 4 — placements
        expand(_OUTDIR + "/4_epang/{batch}/epa_result.jplace", batch=list(_BATCHES.keys())),
        # Step 5a — gappa graft
        expand(_OUTDIR + "/5_gappa/{batch}/graft.newick", batch=list(_BATCHES.keys())),
        # Step 5b — gappa assign
        expand(_OUTDIR + "/5_gappa/{batch}/per_query.tsv", batch=list(_BATCHES.keys())),
        *([
            _OUTDIR + "/6_fasttree/rt_outgroup.faa",
            _OUTDIR + "/6_fasttree/query_with_rt.faa",
            _OUTDIR + "/6_fasttree/my_tree.nwk",
            _OUTDIR + "/6_fasttree/itol_source.txt",
            _OUTDIR + "/6_fasttree/itol_tier.txt",
        ] if _INCLUDE_FASTTREE else []),


# ── Step 0a: extract RdRp_id list for clustered contigs ───────────────────────
# rna_virus_contig_table.tsv has contig_id (col1). tier_summary.tsv has RdRp_id
# (col1) and contig (col3). FAA sequence headers match RdRp_id (e.g.
# AnQing_D_33_0000000012_frame=3), not bare contig_id — so we join on contig to
# get the RdRp_ids and use those as the seqkit grep list.

_TIER_SUMMARY = "result/03_RDRP_identification/4_motif_search/motif_results/tier_summary.tsv"

rule extract_clustered_rdrp_ids:
    input:
        table   = ancient(_RNA_VIRUS_TABLE),
        tier    = ancient(_TIER_SUMMARY)
    output:
        ids = _OUTDIR + "/0b_filtered/clustered_rdrp_ids.txt"
    log:
        err = "log/07b_cluster_RDRP_tree_construction/0b_filtered/extract_ids.err"
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
        mkdir -p $(dirname {output.ids}) $(dirname {log.err})
        python - <<'EOF'
import pandas as pd
table = pd.read_csv("{input.table}", sep="\\t", usecols=["contig_id"])
tier  = pd.read_csv("{input.tier}",  sep="\\t", usecols=["RdRp_id", "contig"])
clustered = set(table["contig_id"])
ids = tier[tier["contig"].isin(clustered)]["RdRp_id"]
ids.to_csv("{output.ids}", index=False, header=False)
print(f"Clustered contigs: {{len(clustered)}}, RdRp_ids written: {{len(ids)}}", flush=True)
EOF
        2> {log.err}
        echo "RdRp_id count: $(wc -l < {output.ids})" >> {log.err}
        """


# ── Step 0b: filter each tier FAA to clustered contigs only ───────────────────
# seqkit grep -f matches sequence names against the RdRp_id list.
# ICTV/ESvirtu FAA files are NOT filtered — they go directly to 0_inputdata.

_FAA_SOURCES = {
    "tier1a_ABCD":            _TIER1A_FAA,
    "tier1b_CABD_depermuted": _TIER1B_FAA,
    "tier2_3motif":           _TIER2_FAA,
}

rule filter_contig_faa:
    input:
        faa = lambda wc: _FAA_SOURCES[wc.tier],
        ids = _OUTDIR + "/0b_filtered/clustered_rdrp_ids.txt"
    output:
        faa = _OUTDIR + "/0b_filtered/contig_{tier}.faa"
    log:
        err = "log/07b_cluster_RDRP_tree_construction/0b_filtered/{tier}.err"
    conda:
        "../envs/mafft_hmmer_seqkit.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p $(dirname {output.faa}) $(dirname {log.err})
        seqkit grep -f {input.ids} {input.faa} > {output.faa} 2> {log.err}
        echo "Sequences retained in {wildcards.tier}: $(grep -c '^>' {output.faa} || echo 0)" >> {log.err}
        """


# ── Step 0: stage filtered contig FAA + optional ICTV/ESvirtu into 0_inputdata ─

rule stage_inputs:
    input:
        f = lambda wc: _STAGED_INPUTS[wc.f],
    output:
        _OUTDIR + "/0_inputdata/{f}",
    log:
        err = "log/07b_cluster_RDRP_tree_construction/0_inputdata/{f}.err",
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log.err})
        cp {input.f} {output} 2> {log.err}
        """


# ── Step 1: hmmbuild profile HMM from the fixed reference MSA ─────────────────

rule hmmbuild_reference:
    input:
        ref_msa = _REF_MSA,
    output:
        hmm = _OUTDIR + "/1_ref_hmm/RCR90.hmm",
    log:
        log = "log/07b_cluster_RDRP_tree_construction/1_ref_hmm/hmmbuild.log",
        err = "log/07b_cluster_RDRP_tree_construction/1_ref_hmm/hmmbuild.err",
    conda:
        "../envs/hmmer.yaml"
    threads: config["epa_ng"]["hmmbuild"]["threads"]
    resources:
        mem_mb_per_cpu  = config["epa_ng"]["hmmbuild"]["memory"],
        runtime         = config["epa_ng"]["hmmbuild"]["runtime"],
        cpus_per_task   = config["epa_ng"]["hmmbuild"]["threads"],
        slurm_partition = config["epa_ng"]["hmmbuild"]["partition"],
        slurm_account   = config["epa_ng"]["hmmbuild"]["account"],
    shell:
        """
        mkdir -p $(dirname {output.hmm}) $(dirname {log.log})
        hmmbuild --amino --symfrac 0.0 --cpu {threads} \
            {output.hmm} {input.ref_msa} \
            > {log.log} 2> {log.err}
        """


# ── Step 2: hmmalign --trim queries into reference column coordinates ────────

rule hmmalign_query:
    input:
        hmm   = _OUTDIR + "/1_ref_hmm/RCR90.hmm",
        query = lambda wc: _BATCHES[wc.batch],
    output:
        sto = _OUTDIR + "/2_hmmalign/{batch}/aligned.sto",
    params:
        merged_query = _OUTDIR + "/2_hmmalign/{batch}/query_merged.faa",
    log:
        log = "log/07b_cluster_RDRP_tree_construction/2_hmmalign/{batch}.log",
        err = "log/07b_cluster_RDRP_tree_construction/2_hmmalign/{batch}.err",
    conda:
        "../envs/hmmer.yaml"
    threads: config["epa_ng"]["hmmalign"]["threads"]
    resources:
        mem_mb_per_cpu  = config["epa_ng"]["hmmalign"]["memory"],
        runtime         = config["epa_ng"]["hmmalign"]["runtime"],
        cpus_per_task   = config["epa_ng"]["hmmalign"]["threads"],
        slurm_partition = config["epa_ng"]["hmmalign"]["partition"],
        slurm_account   = config["epa_ng"]["hmmalign"]["account"],
    shell:
        """
        mkdir -p $(dirname {output.sto}) $(dirname {log.log})
        cat {input.query} > {params.merged_query}
        hmmalign --trim --amino --outformat Stockholm \
            -o {output.sto} \
            {input.hmm} {params.merged_query} \
            > {log.log} 2> {log.err}
        """


# ── Step 3: split aligned Stockholm into reference-MSA and query-MSA FASTA ───

rule split_ref_query:
    input:
        sto     = _OUTDIR + "/2_hmmalign/{batch}/aligned.sto",
        ref_msa = _REF_MSA,
    output:
        ref_fasta   = _OUTDIR + "/3_split/{batch}/reference.fasta",
        query_fasta = _OUTDIR + "/3_split/{batch}/query.fasta",
    log:
        err = "log/07b_cluster_RDRP_tree_construction/3_split/{batch}.err",
    conda:
        "../envs/mafft_hmmer_seqkit.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p $(dirname {output.ref_fasta})
        esl-reformat afa {input.ref_msa} > {output.ref_fasta} 2> {log.err}
        esl-alimask --rf-is-mask {input.sto} 2>> {log.err} \
            | esl-reformat afa - > {output.query_fasta} 2>> {log.err}
        """


# ── Step 4: EPA-ng placement onto the FIXED reference tree ───────────────────

rule epa_ng_place:
    input:
        tree    = _REF_TREE,
        ref_msa = _OUTDIR + "/3_split/{batch}/reference.fasta",
        query   = _OUTDIR + "/3_split/{batch}/query.fasta",
    output:
        jplace = _OUTDIR + "/4_epang/{batch}/epa_result.jplace",
    params:
        outdir = _OUTDIR + "/4_epang/{batch}",
        model  = _MODEL,
    log:
        log = "log/07b_cluster_RDRP_tree_construction/4_epang/{batch}.log",
        err = "log/07b_cluster_RDRP_tree_construction/4_epang/{batch}.err",
    conda:
        "../envs/epang_gappa.yaml"
    threads: config["epa_ng"]["epang"]["threads"]
    resources:
        mem_mb_per_cpu  = config["epa_ng"]["epang"]["memory"],
        runtime         = config["epa_ng"]["epang"]["runtime"],
        cpus_per_task   = config["epa_ng"]["epang"]["threads"],
        slurm_partition = config["epa_ng"]["epang"]["partition"],
        slurm_account   = config["epa_ng"]["epang"]["account"],
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log.log})
        epa-ng \
            --tree {input.tree} \
            --ref-msa {input.ref_msa} \
            --query {input.query} \
            --model {params.model} \
            --outdir {params.outdir} \
            --redo \
            -T {threads} \
            > {log.log} 2> {log.err}
        """


# ── Step 5a: gappa graft ──────────────────────────────────────────────────────

rule gappa_graft:
    input:
        jplace = _OUTDIR + "/4_epang/{batch}/epa_result.jplace",
    output:
        newick = _OUTDIR + "/5_gappa/{batch}/graft.newick",
    params:
        outdir = _OUTDIR + "/5_gappa/{batch}",
    log:
        log = "log/07b_cluster_RDRP_tree_construction/5_gappa/{batch}.graft.log",
        err = "log/07b_cluster_RDRP_tree_construction/5_gappa/{batch}.graft.err",
    conda:
        "../envs/epang_gappa.yaml"
    threads: config["epa_ng"]["gappa"]["threads"]
    resources:
        mem_mb_per_cpu  = config["epa_ng"]["gappa"]["memory"],
        runtime         = config["epa_ng"]["gappa"]["runtime"],
        cpus_per_task   = config["epa_ng"]["gappa"]["threads"],
        slurm_partition = config["epa_ng"]["gappa"]["partition"],
        slurm_account   = config["epa_ng"]["gappa"]["account"],
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log.log})
        gappa examine graft \
            --jplace-path {input.jplace} \
            --out-dir {params.outdir} \
            --allow-file-overwriting \
            > {log.log} 2> {log.err}
        mv {params.outdir}/*.newick {output.newick} 2>/dev/null || true
        """


# ══════════════════════════════════════════════════════════════════════════════
# Taxonomy preparation + gappa assign
# ══════════════════════════════════════════════════════════════════════════════

rule prepare_taxon_file:
    input:
        tsv = _RVMT_INFO_TSV,
    output:
        tsv = _TAXON_FILE,
    log:
        err = "log/07b_cluster_RDRP_tree_construction/0_taxon/prepare_taxon_file.err",
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log.err})
        awk -F'\\t' 'NR>1 && $3 != "" && $3 !~ /^rt\\./ && $12 ~ /^Lvl 0/ {{
            print $3 "\\t" $16 ";" $17 ";" $18 ";" $19 ";" $20
        }}' {input.tsv} > {output.tsv} 2> {log.err}
        echo "Taxon entries written: $(wc -l < {output.tsv})" >> {log.err}
        """


rule gappa_assign:
    input:
        jplace     = _OUTDIR + "/4_epang/{batch}/epa_result.jplace",
        taxon_file = _TAXON_FILE,
    output:
        tsv = _OUTDIR + "/5_gappa/{batch}/per_query.tsv",
    params:
        outdir = _OUTDIR + "/5_gappa/{batch}",
    log:
        log = "log/07b_cluster_RDRP_tree_construction/5_gappa/{batch}.assign.log",
        err = "log/07b_cluster_RDRP_tree_construction/5_gappa/{batch}.assign.err",
    conda:
        "../envs/epang_gappa.yaml"
    threads: config["epa_ng"]["gappa"]["threads"]
    resources:
        mem_mb_per_cpu  = config["epa_ng"]["gappa"]["memory"],
        runtime         = config["epa_ng"]["gappa"]["runtime"],
        cpus_per_task   = config["epa_ng"]["gappa"]["threads"],
        slurm_partition = config["epa_ng"]["gappa"]["partition"],
        slurm_account   = config["epa_ng"]["gappa"]["account"],
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log.log})
        gappa examine assign \
            --jplace-path {input.jplace} \
            --taxon-file {input.taxon_file} \
            --out-dir {params.outdir} \
            --allow-file-overwriting \
            --per-query-results \
            > {log.log} 2> {log.err}
        mv {params.outdir}/per_query.tsv {output.tsv} 2>/dev/null || true
        """


if _INCLUDE_FASTTREE:
    rule extract_rt_outgroup:
        input:
            ref_msa = _REF_MSA,
        output:
            faa = _OUTDIR + "/6_fasttree/rt_outgroup.faa",
        log:
            err = "log/07b_cluster_RDRP_tree_construction/6_fasttree/extract_rt_outgroup.err",
        conda:
            "../envs/mafft_hmmer_seqkit.yaml"
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.faa}) $(dirname {log.err})
            esl-reformat afa {input.ref_msa} 2>> {log.err} \
                | seqkit grep --use-regexp --pattern '^rt\\.' --threads {threads} \
                > {output.faa} 2>> {log.err}
            echo "RT outgroup count: $(grep -c '^>' {output.faa})" >> {log.err}
            """

    rule concat_query_rt:
        input:
            query = _OUTDIR + "/3_split/tier1ab2/query.fasta",
            rt    = _OUTDIR + "/6_fasttree/rt_outgroup.faa",
        output:
            faa = _OUTDIR + "/6_fasttree/query_with_rt.faa",
        log:
            err = "log/07b_cluster_RDRP_tree_construction/6_fasttree/concat_query_rt.err",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.faa}) $(dirname {log.err})
            cat {input.query} {input.rt} > {output.faa} 2> {log.err}
            echo "Total sequences: $(grep -c '^>' {output.faa})" >> {log.err}
            """

    rule fasttree_build:
        input:
            faa = _OUTDIR + "/6_fasttree/query_with_rt.faa",
        output:
            nwk = _OUTDIR + "/6_fasttree/my_tree.nwk",
        log:
            log = "log/07b_cluster_RDRP_tree_construction/6_fasttree/fasttree.log",
            err = "log/07b_cluster_RDRP_tree_construction/6_fasttree/fasttree.err",
        conda:
            "../envs/fasttree.yaml"
        threads: config["epa_ng"]["fasttree"]["threads"]
        resources:
            mem_mb_per_cpu  = config["epa_ng"]["fasttree"]["memory"],
            runtime         = config["epa_ng"]["fasttree"]["runtime"],
            cpus_per_task   = config["epa_ng"]["fasttree"]["threads"],
            slurm_partition = config["epa_ng"]["fasttree"]["partition"],
            slurm_account   = config["epa_ng"]["fasttree"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.nwk}) $(dirname {log.log})
            export OMP_NUM_THREADS={threads}
            FastTree -wag -gamma \
                {input.faa} \
                > {output.nwk} \
                2> {log.log}
            """

    rule itol_colorstrip:
        input:
            faa          = _OUTDIR + "/6_fasttree/query_with_rt.faa",
            tier_summary = ancient(_MOTIF_RESULTS + "/tier_summary.tsv"),
        output:
            itol_source = _OUTDIR + "/6_fasttree/itol_source.txt",
            itol_tier   = _OUTDIR + "/6_fasttree/itol_tier.txt",
        params:
            script = "scripts/07_phylogenetic_tree/itol_colorstrip.py",
        log:
            err = "log/07b_cluster_RDRP_tree_construction/6_fasttree/itol_colorstrip.err",
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
            python {params.script} \
                --faa          {input.faa} \
                --tier-summary {input.tier_summary} \
                --output       {output.itol_source} \
                --output-tier  {output.itol_tier} \
                2> {log.err}
            """
