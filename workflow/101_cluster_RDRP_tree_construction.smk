# RdRP phylogenetic PLACEMENT pipeline — filtered contigs only (after workflow 99)
#
# Standalone pipeline; does NOT depend on 07_RDRP_tree_construction.smk.
# Inputs come from workflow 03 (final confirmed FAA files) and workflow 99
# (rna_virus_contig_table.tsv).
#
# Differences from 07:
#   - Only contigs present in rna_virus_contig_table.tsv are used as queries.
#     The nr_filtered/RdRpCATCH.faa / LucaProt.faa files are NR-confirmed in
#     Step 0b to retain only those contig IDs.
#   - ICTV / ESvirtu reference candidates are still included when enabled in
#     config, because they are not assembly contigs and are not filtered.
#   - Outputs go to result/101_cluster_RDRP_tree_construction/
#   - Log files go to log/101_cluster_RDRP_tree_construction/
#
# Run AFTER 99_tables.smk has produced rna_virus_contig_table.tsv.

import os

configfile: "config/config.yaml"

_FINAL_DIR   = "result/03_RDRP_identification/10_final"
_RC_FAA      = _FINAL_DIR + "/nr_filtered/RdRpCATCH.faa"
_LP_FAA      = _FINAL_DIR + "/nr_filtered/LucaProt.faa"

_ICTV_PROTEIN_DIR  = "result/03_RDRP_identification/ICTV/8_RdRp_protein"
_ICTV_RC_FAA       = _ICTV_PROTEIN_DIR + "/RdRpCATCH.faa"   # full ICTV proteins identified by RdRpCATCH
_ICTV_LP_FAA       = _ICTV_PROTEIN_DIR + "/LucaProt.faa"    # full ICTV proteins identified by LucaProt

_ESV_PROTEIN_DIR   = "result/03b_esvirtu_rdrp_identification/10_proteins"
_ESV_RC_FAA        = _ESV_PROTEIN_DIR + "/RdRpCATCH.faa"
_ESV_LP_FAA        = _ESV_PROTEIN_DIR + "/LucaProt.faa"

# Clustered contig table produced by workflow 99
_RNA_VIRUS_TABLE   = "result/99_tables/4_contig_tag_abundance/rna_virus_contig_table.tsv"
_RDRP_MERGED_TSV   = _FINAL_DIR + "/final_rdrp_merged.tsv"
_ESVIRTU_INFO_TSV  = "result/99_tables/2_esvirtu_table/esvirtu_info_table.tsv"
_SPEARMAN_DIR      = "result/99_tables/6_Recalculated_correlation/1_RNA_virus_Recalculated_correlation_TPM/2_spearman_analysis"
_ICTV_VMR_XLSX     = config["ICTV"]["vmr_xlsx"]
_ICTV_VMR_SHEET    = config["ICTV"]["vmr_sheet"]
_ITOL_MIN_SEED_NONZERO   = config["epa_ng"].get("itol", {}).get("min_seed_nonzero",   0)
_ITOL_MIN_TARGET_NONZERO = config["epa_ng"].get("itol", {}).get("min_target_nonzero", 0)

_OUTDIR      = "result/101_cluster_RDRP_tree_construction"
_REF_MSA     = config["epa_ng"]["ref_msa"]
_REF_TREE    = config["epa_ng"]["ref_tree"]
_MODEL       = config["epa_ng"]["model"]
_RVMT_INFO_TSV = config["epa_ng"]["rvmt_info_tsv"]
_TAXON_FILE  = _OUTDIR + "/0_taxon/RCR90_taxonomy.tsv"

_INCLUDE_ICTV           = config["epa_ng"].get("include_ictv", False)
_INCLUDE_ESVIRTU        = config["epa_ng"].get("include_esvirtu", False)
_FASTTREE_BOOTSTRAP = config["epa_ng"]["fasttree"].get("bootstrap", 0)
_FASTTREE_BOOT_FLAG = f"-boot {_FASTTREE_BOOTSTRAP}" if _FASTTREE_BOOTSTRAP > 0 else ""

_INCLUDE_FASTTREE_FILTERED      = config["epa_ng"].get("include_fasttree", False)
_INCLUDE_FASTTREE_FILTERED_ICTV = config["epa_ng"].get("include_fasttree_ictv", False)
_INCLUDE_FASTTREE_UNFILTERED      = config["epa_ng"].get("include_fasttree_full", False)
_INCLUDE_FASTTREE_UNFILTERED_ICTV = config["epa_ng"].get("include_fasttree_full_ictv", False)
_INCLUDE_EPANG              = config["epa_ng"].get("include_epang", True)
_INCLUDE_EPANG_UNFILTERED   = config["epa_ng"].get("include_epang_unfiltered", False)

# Use epa_ng-specific small_job so partition stays on "participation"
# without affecting other workflows that use the global small_job block.
_SJ = config["epa_ng"]["small_job"]

_CP_TABLE_DIR = "result/99_tables/01_cp_table/101_cluster_RDRP_tree_construction"

# Filtered contig FAA files (0b_filtered); ICTV/ESvirtu staged directly
_STAGED_INPUTS = {
    "contig_RdRpCATCH.faa":       _OUTDIR + "/0b_filtered/contig_RdRpCATCH.faa",
    "contig_LucaProt.faa":         _OUTDIR + "/0b_filtered/contig_LucaProt.faa",
    "unfiltered_RdRpCATCH.faa":   _RC_FAA,
    "unfiltered_LucaProt.faa":    _LP_FAA,
}
if _INCLUDE_ICTV or _INCLUDE_FASTTREE_FILTERED_ICTV or _INCLUDE_FASTTREE_UNFILTERED_ICTV:
    _STAGED_INPUTS.update({
        "ictv_RdRpCATCH.faa": _ICTV_RC_FAA,
        "ictv_LucaProt.faa":  _ICTV_LP_FAA,
    })
if _INCLUDE_ESVIRTU:
    _STAGED_INPUTS.update({
        "esvirtu_RdRpCATCH.faa": _ESV_RC_FAA,
        "esvirtu_LucaProt.faa":  _ESV_LP_FAA,
    })

# Batch definitions — each batch goes through hmmalign_query + split_ref_query
#   filtered:          filtered contigs [+ ESvirtu]               → EPA-ng + filtered tree
#   filtered_ictv:     filtered contigs [+ ESvirtu] + ICTV        → filtered+ICTV tree
#   unfiltered:        all 4010+142 proteins + ESvirtu             → unfiltered tree
#   unfiltered_ictv:   all 4010+142 proteins + ESvirtu + ICTV     → unfiltered+ICTV tree
_CONTIG_FAAS = [_OUTDIR + "/0_inputdata/contig_RdRpCATCH.faa",
                _OUTDIR + "/0_inputdata/contig_LucaProt.faa"]
_ESVIRTU_FAAS = ([_OUTDIR + "/0_inputdata/esvirtu_RdRpCATCH.faa",
                  _OUTDIR + "/0_inputdata/esvirtu_LucaProt.faa"] if _INCLUDE_ESVIRTU else [])
_NEED_ICTV = _INCLUDE_ICTV or _INCLUDE_FASTTREE_FILTERED_ICTV or _INCLUDE_FASTTREE_UNFILTERED_ICTV
_ICTV_FAAS = ([_OUTDIR + "/0_inputdata/ictv_RdRpCATCH.faa",
               _OUTDIR + "/0_inputdata/ictv_LucaProt.faa"] if _NEED_ICTV else [])
_UNFILTERED_FAAS = [_OUTDIR + "/0_inputdata/unfiltered_RdRpCATCH.faa",
                    _OUTDIR + "/0_inputdata/unfiltered_LucaProt.faa"]

_BATCHES = {
    "filtered": _CONTIG_FAAS + _ESVIRTU_FAAS,
}
if _INCLUDE_FASTTREE_FILTERED_ICTV:
    _BATCHES["filtered_ictv"] = _CONTIG_FAAS + _ESVIRTU_FAAS + _ICTV_FAAS
if _INCLUDE_FASTTREE_UNFILTERED:
    _BATCHES["unfiltered"] = _UNFILTERED_FAAS + _ESVIRTU_FAAS
if _INCLUDE_EPANG_UNFILTERED:
    # separate batch for EPA-ng: unfiltered contigs only, no ESvirtu
    _BATCHES["unfiltered_epang"] = _UNFILTERED_FAAS
if _INCLUDE_FASTTREE_UNFILTERED_ICTV:
    _BATCHES["unfiltered_ictv"] = _UNFILTERED_FAAS + _ESVIRTU_FAAS + _ICTV_FAAS

# EPA-ng runs only on batches explicitly meant for placement (not FastTree-only batches)
_EPANG_BATCHES = {}
if _INCLUDE_EPANG:
    _EPANG_BATCHES["filtered"] = _BATCHES["filtered"]
if _INCLUDE_EPANG_UNFILTERED:
    _EPANG_BATCHES["unfiltered_epang"] = _BATCHES["unfiltered_epang"]


rule all:
    input:
        _TAXON_FILE,
        _OUTDIR + "/0b_filtered/clustered_rdrp_ids.txt",
        expand(_OUTDIR + "/0b_filtered/contig_{cat}.faa", cat=["RdRpCATCH", "LucaProt"]),
        expand(_OUTDIR + "/0_inputdata/{f}", f=list(_STAGED_INPUTS.keys())),
        _OUTDIR + "/1_ref_hmm/RCR90.hmm",
        expand(_OUTDIR + "/2_hmmalign/{batch}/aligned.sto", batch=list(_BATCHES.keys())),
        expand(_OUTDIR + "/3_split/{batch}/query.fasta",     batch=list(_BATCHES.keys())),
        expand(_OUTDIR + "/3_split/{batch}/reference.fasta", batch=list(_BATCHES.keys())),
        *expand(_OUTDIR + "/4_epang/{batch}/epa_result.jplace",          batch=list(_EPANG_BATCHES.keys())),
        *expand(_OUTDIR + "/5_gappa/{batch}/graft.newick",               batch=list(_EPANG_BATCHES.keys())),
        *expand(_OUTDIR + "/5_gappa/{batch}/per_query.tsv",              batch=list(_EPANG_BATCHES.keys())),
        *expand(_OUTDIR + "/5_gappa/{batch}/phylum_annotation.tsv",      batch=list(_EPANG_BATCHES.keys())),
        *expand(_OUTDIR + "/5_gappa/{batch}/full_taxonomy_annotation.tsv", batch=list(_EPANG_BATCHES.keys())),
        *([
            _OUTDIR + "/6_fasttree/rt_outgroup.faa",
            _OUTDIR + "/6_fasttree/filtered/query_with_rt.faa",
            _OUTDIR + "/6_fasttree/filtered/my_tree.nwk",
            _OUTDIR + "/6_fasttree/filtered/my_tree.renamed.nwk",
            _OUTDIR + "/6_fasttree/filtered/itol/colorstrip/itol_category.txt",
            _OUTDIR + "/6_fasttree/filtered/itol/colorstrip/itol_phylum_gappa.txt",
            *expand(_OUTDIR + "/6_fasttree/filtered/itol/connection_by_rank/itol_connection_r{t}.txt",
                    t=["0.6", "0.7", "0.8", "0.9"]),
            *expand(_OUTDIR + "/6_fasttree/filtered/itol/connection_by_virus/itol_connection_r{t}.txt",
                    t=["0.6", "0.7", "0.8", "0.9"]),
        ] if _INCLUDE_FASTTREE_FILTERED else []),
        *([
            _OUTDIR + "/6_fasttree/filtered_ictv/query_with_rt.faa",
            _OUTDIR + "/6_fasttree/filtered_ictv/my_tree.nwk",
            _OUTDIR + "/6_fasttree/filtered_ictv/my_tree.renamed.nwk",
            _OUTDIR + "/6_fasttree/filtered_ictv/contig_classification/contig_classification.tsv",
            *expand(_OUTDIR + "/6_fasttree/filtered_ictv/itol/colorstrip/ictv/itol_ictv_{rank}.txt",
                    rank=["phylum", "class", "order", "family", "genus"]),
        ] if _INCLUDE_FASTTREE_FILTERED_ICTV else []),
        *([
            _OUTDIR + "/6_fasttree/unfiltered/query_with_rt.faa",
            _OUTDIR + "/6_fasttree/unfiltered/my_tree.nwk",
            _OUTDIR + "/6_fasttree/unfiltered/my_tree.renamed.nwk",
        ] if _INCLUDE_FASTTREE_UNFILTERED else []),
        *([
            _OUTDIR + "/6_fasttree/unfiltered_ictv/query_with_rt.faa",
            _OUTDIR + "/6_fasttree/unfiltered_ictv/my_tree.nwk",
            _OUTDIR + "/6_fasttree/unfiltered_ictv/my_tree.renamed.nwk",
        ] if _INCLUDE_FASTTREE_UNFILTERED_ICTV else []),
        # ── Copy taxonomy annotations to 99_tables ────────────────────────────
        *expand(_CP_TABLE_DIR + "/5_gappa/{batch}/phylum_annotation.tsv",       batch=list(_EPANG_BATCHES.keys())),
        *expand(_CP_TABLE_DIR + "/5_gappa/{batch}/full_taxonomy_annotation.tsv", batch=list(_EPANG_BATCHES.keys())),
        *([_CP_TABLE_DIR + "/6_fasttree/filtered_ictv/contig_classification.tsv"]
          if _INCLUDE_FASTTREE_FILTERED_ICTV else []),


# ── Step 0a: extract protein IDs for contigs in rna_virus_contig_table ────────
# final_rdrp_merged.tsv has contig (col2) and source (RdRpCATCH/LucaProt).
# FAA headers match the protein_id format (e.g. GouBa_N_1_0000000004_frame=-1_RdRp_27-1377
# for RdRpCATCH, ORF1_GouBa_N_1_0000000399:1:1680 for LucaProt).
# We extract all protein IDs whose contig appears in rna_virus_contig_table.tsv.

rule extract_clustered_rdrp_ids:
    input:
        table  = ancient(_RNA_VIRUS_TABLE),
        merged = ancient(_RDRP_MERGED_TSV),
        rc_faa = ancient(_RC_FAA),
        lp_faa = ancient(_LP_FAA),
    output:
        ids = _OUTDIR + "/0b_filtered/clustered_rdrp_ids.txt",
    log:
        err = "log/101_cluster_RDRP_tree_construction/0b_filtered/extract_ids.err",
    conda:
        "../envs/python.yaml"
    threads: _SJ["threads"]
    resources:
        mem_mb_per_cpu  = _SJ["memory"],
        runtime         = _SJ["runtime"],
        cpus_per_task   = _SJ["threads"],
        slurm_partition = _SJ["partition"],
        slurm_account   = _SJ["account"],
    script:
        "../scripts/07_phylogenetic_tree/extract_clustered_rdrp_ids.py"


# ── Step 0b: filter each FAA to filtered contigs only ────────────────────────

_FAA_SOURCES = {
    "RdRpCATCH": _RC_FAA,
    "LucaProt":  _LP_FAA,
}

rule filter_contig_faa:
    input:
        faa = lambda wc: _FAA_SOURCES[wc.cat],
        ids = _OUTDIR + "/0b_filtered/clustered_rdrp_ids.txt"
    output:
        faa = _OUTDIR + "/0b_filtered/contig_{cat}.faa"
    wildcard_constraints:
        cat = "RdRpCATCH|LucaProt"
    log:
        err = "log/101_cluster_RDRP_tree_construction/0b_filtered/{cat}.err"
    conda:
        "../envs/mafft_hmmer_seqkit.yaml"
    threads: _SJ["threads"]
    resources:
        mem_mb_per_cpu  = _SJ["memory"],
        runtime         = _SJ["runtime"],
        cpus_per_task   = _SJ["threads"],
        slurm_partition = _SJ["partition"],
        slurm_account   = _SJ["account"],
    shell:
        """
        mkdir -p $(dirname {output.faa}) $(dirname {log.err})
        seqkit grep -f {input.ids} {input.faa} > {output.faa} 2> {log.err}
        echo "Sequences retained ({wildcards.cat}): $(grep -c '^>' {output.faa} || echo 0)" >> {log.err}
        """


# ── Step 0: stage filtered contig FAA + optional ICTV/ESvirtu into 0_inputdata ─

rule stage_inputs:
    input:
        f = lambda wc: _STAGED_INPUTS[wc.f],
    output:
        _OUTDIR + "/0_inputdata/{f}",
    wildcard_constraints:
        f = "(?!ictv_).+",
    log:
        err = "log/101_cluster_RDRP_tree_construction/0_inputdata/{f}.err",
    threads: _SJ["threads"]
    resources:
        mem_mb_per_cpu  = _SJ["memory"],
        runtime         = _SJ["runtime"],
        cpus_per_task   = _SJ["threads"],
        slurm_partition = _SJ["partition"],
        slurm_account   = _SJ["account"],
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log.err})
        cp {input.f} {output} 2> {log.err}
        """

rule stage_ictv_inputs:
    input:
        f = lambda wc: _STAGED_INPUTS[wc.f],
    output:
        _OUTDIR + "/0_inputdata/{f}",
    wildcard_constraints:
        f = "ictv_.+",
    log:
        err = "log/101_cluster_RDRP_tree_construction/0_inputdata/{f}.err",
    threads: _SJ["threads"]
    resources:
        mem_mb_per_cpu  = _SJ["memory"],
        runtime         = _SJ["runtime"],
        cpus_per_task   = _SJ["threads"],
        slurm_partition = _SJ["partition"],
        slurm_account   = _SJ["account"],
    shell:
        """
        mkdir -p $(dirname {output}) $(dirname {log.err})
        sed 's/^>/\\&ICTV_/' {input.f} > {output} 2> {log.err}
        """


# ── Step 1: hmmbuild profile HMM from the fixed reference MSA ─────────────────

rule hmmbuild_reference:
    input:
        ref_msa = _REF_MSA,
    output:
        hmm = _OUTDIR + "/1_ref_hmm/RCR90.hmm",
    log:
        log = "log/101_cluster_RDRP_tree_construction/1_ref_hmm/hmmbuild.log",
        err = "log/101_cluster_RDRP_tree_construction/1_ref_hmm/hmmbuild.err",
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


# ── Step 2: hmmalign --trim queries into reference column coordinates ──────────

rule hmmalign_query:
    input:
        hmm   = _OUTDIR + "/1_ref_hmm/RCR90.hmm",
        query = lambda wc: _BATCHES[wc.batch],
    output:
        sto = _OUTDIR + "/2_hmmalign/{batch}/aligned.sto",
    params:
        merged_query = _OUTDIR + "/2_hmmalign/{batch}/query_merged.faa",
    log:
        log = "log/101_cluster_RDRP_tree_construction/2_hmmalign/{batch}.log",
        err = "log/101_cluster_RDRP_tree_construction/2_hmmalign/{batch}.err",
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


# ── Step 3: split aligned Stockholm into reference-MSA and query-MSA FASTA ────

rule split_ref_query:
    input:
        sto     = _OUTDIR + "/2_hmmalign/{batch}/aligned.sto",
        ref_msa = _REF_MSA,
    output:
        ref_fasta   = _OUTDIR + "/3_split/{batch}/reference.fasta",
        query_fasta = _OUTDIR + "/3_split/{batch}/query.fasta",
    log:
        err = "log/101_cluster_RDRP_tree_construction/3_split/{batch}.err",
    conda:
        "../envs/mafft_hmmer_seqkit.yaml"
    threads: _SJ["threads"]
    resources:
        mem_mb_per_cpu  = _SJ["memory"],
        runtime         = _SJ["runtime"],
        cpus_per_task   = _SJ["threads"],
        slurm_partition = _SJ["partition"],
        slurm_account   = _SJ["account"],
    shell:
        """
        mkdir -p $(dirname {output.ref_fasta})
        esl-reformat afa {input.ref_msa} > {output.ref_fasta} 2> {log.err}
        esl-alimask --rf-is-mask {input.sto} 2>> {log.err} \
            | esl-reformat afa - > {output.query_fasta} 2>> {log.err}
        """


# ── Step 4-5: EPA-ng placement + gappa (conditional) ─────────────────────────

if _EPANG_BATCHES:
    rule stage_epang_inputs:
        wildcard_constraints:
            batch = "|".join(_EPANG_BATCHES.keys()),
        input:
            faa = lambda wc: _EPANG_BATCHES[wc.batch],
        output:
            done = _OUTDIR + "/4_epang/0_input_file/{batch}/.done",
        params:
            outdir = _OUTDIR + "/4_epang/0_input_file/{batch}",
        log:
            err = "log/101_cluster_RDRP_tree_construction/4_epang/0_input_file/{batch}.err",
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p {params.outdir} $(dirname {log.err})
            for f in {input.faa}; do
                cp "$f" {params.outdir}/$(basename "$f") 2>> {log.err}
            done
            touch {output.done}
            """

    rule epa_ng_place:
        wildcard_constraints:
            batch = "|".join(_EPANG_BATCHES.keys()),
        input:
            tree    = _REF_TREE,
            ref_msa = _OUTDIR + "/3_split/{batch}/reference.fasta",
            query   = _OUTDIR + "/3_split/{batch}/query.fasta",
            staged  = _OUTDIR + "/4_epang/0_input_file/{batch}/.done",
        output:
            jplace = _OUTDIR + "/4_epang/{batch}/epa_result.jplace",
        params:
            outdir = _OUTDIR + "/4_epang/{batch}",
            model  = _MODEL,
        log:
            log = "log/101_cluster_RDRP_tree_construction/4_epang/{batch}.log",
            err = "log/101_cluster_RDRP_tree_construction/4_epang/{batch}.err",
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

    rule gappa_graft:
        wildcard_constraints:
            batch = "|".join(_EPANG_BATCHES.keys()),
        input:
            jplace = _OUTDIR + "/4_epang/{batch}/epa_result.jplace",
        output:
            newick = _OUTDIR + "/5_gappa/{batch}/graft.newick",
        params:
            outdir = _OUTDIR + "/5_gappa/{batch}",
        log:
            log = "log/101_cluster_RDRP_tree_construction/5_gappa/{batch}.graft.log",
            err = "log/101_cluster_RDRP_tree_construction/5_gappa/{batch}.graft.err",
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

    rule gappa_assign:
        wildcard_constraints:
            batch = "|".join(_EPANG_BATCHES.keys()),
        input:
            jplace     = _OUTDIR + "/4_epang/{batch}/epa_result.jplace",
            taxon_file = _TAXON_FILE,
        output:
            tsv = _OUTDIR + "/5_gappa/{batch}/per_query.tsv",
        params:
            outdir = _OUTDIR + "/5_gappa/{batch}",
        log:
            log = "log/101_cluster_RDRP_tree_construction/5_gappa/{batch}.assign.log",
            err = "log/101_cluster_RDRP_tree_construction/5_gappa/{batch}.assign.err",
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

    rule parse_gappa_taxonomy:
        wildcard_constraints:
            batch = "|".join(_EPANG_BATCHES.keys()),
        input:
            tsv    = _OUTDIR + "/5_gappa/{batch}/per_query.tsv",
            script = "scripts/07_phylogenetic_tree/parse_gappa_per_query.py",
        output:
            phylum = _OUTDIR + "/5_gappa/{batch}/phylum_annotation.tsv",
            full   = _OUTDIR + "/5_gappa/{batch}/full_taxonomy_annotation.tsv",
        params:
            outdir    = _OUTDIR + "/5_gappa/{batch}",
            threshold = config["epa_ng"].get("gappa_afract_threshold", 0.66),
        log:
            err = "log/101_cluster_RDRP_tree_construction/5_gappa/{batch}.parse_taxonomy.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            python {input.script} \
                --input     {input.tsv} \
                --outdir    {params.outdir} \
                --threshold {params.threshold} \
                2> {log.err}
            """


# ── Taxonomy preparation (always runs — needed by gappa assign) ───────────────

rule prepare_taxon_file:
    input:
        tsv = _RVMT_INFO_TSV,
    output:
        tsv = _TAXON_FILE,
    log:
        err = "log/101_cluster_RDRP_tree_construction/0_taxon/prepare_taxon_file.err",
    threads: _SJ["threads"]
    resources:
        mem_mb_per_cpu  = _SJ["memory"],
        runtime         = _SJ["runtime"],
        cpus_per_task   = _SJ["threads"],
        slurm_partition = _SJ["partition"],
        slurm_account   = _SJ["account"],
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log.err})
        awk -F'\\t' 'NR>1 && $3 != "" && $3 !~ /^rt\\./ && $12 ~ /^Lvl 0/ {{
            print $3 "\\t" $16 ";" $17 ";" $18 ";" $19 ";" $20
        }}' {input.tsv} > {output.tsv} 2> {log.err}
        echo "Taxon entries written: $(wc -l < {output.tsv})" >> {log.err}
        """


if _INCLUDE_FASTTREE_FILTERED:
    rule extract_rt_outgroup:
        input:
            ref_msa = _REF_MSA,
        output:
            faa = _OUTDIR + "/6_fasttree/rt_outgroup.faa",  # shared by all trees
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/extract_rt_outgroup.err",
        conda:
            "../envs/mafft_hmmer_seqkit.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output.faa}) $(dirname {log.err})
            esl-reformat afa {input.ref_msa} 2>> {log.err} \
                | seqkit grep --use-regexp --pattern '^rt\\.' --threads {threads} \
                > {output.faa} 2>> {log.err}
            echo "RT outgroup count: $(grep -c '^>' {output.faa})" >> {log.err}
            """

    rule concat_query_rt_filtered:
        input:
            query = _OUTDIR + "/3_split/filtered/query.fasta",
            rt    = _OUTDIR + "/6_fasttree/rt_outgroup.faa",
        output:
            faa = _OUTDIR + "/6_fasttree/filtered/query_with_rt.faa",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered/concat_query_rt.err",
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output.faa}) $(dirname {log.err})
            cat {input.query} {input.rt} > {output.faa} 2> {log.err}
            echo "Total sequences (filtered): $(grep -c '^>' {output.faa})" >> {log.err}
            """

    rule fasttree_build_filtered:
        input:
            faa = _OUTDIR + "/6_fasttree/filtered/query_with_rt.faa",
        output:
            nwk = _OUTDIR + "/6_fasttree/filtered/my_tree.nwk",
        log:
            log = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered/fasttree.log",
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered/fasttree.err",
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
            FastTree -wag -gamma {_FASTTREE_BOOT_FLAG} \
                {input.faa} \
                > {output.nwk} \
                2> {log.log}
            """

    rule rename_tree_tips_filtered:
        input:
            nwk = _OUTDIR + "/6_fasttree/filtered/my_tree.nwk",
        output:
            nwk    = _OUTDIR + "/6_fasttree/filtered/my_tree.renamed.nwk",
            lookup = _OUTDIR + "/6_fasttree/filtered/my_tree.renamed_lookup.tsv",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered/rename.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {log.err})
            python scripts/07_phylogenetic_tree/rename_tree_tips.py \
                --input  {input.nwk} \
                --output {output.nwk} \
                --lookup {output.lookup} \
                2> {log.err}
            """

    rule itol_colorstrip:
        input:
            faa    = _OUTDIR + "/6_fasttree/filtered/query_with_rt.faa",
            lookup = _OUTDIR + "/6_fasttree/filtered/my_tree.renamed_lookup.tsv",
            merged = ancient(_RDRP_MERGED_TSV),
            script = "scripts/07_phylogenetic_tree/itol_colorstrip.py",
        output:
            itol = _OUTDIR + "/6_fasttree/filtered/itol/colorstrip/itol_category.txt",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered/itol_colorstrip.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output.itol})
            python {input.script} \
                --faa         {input.faa} \
                --lookup      {input.lookup} \
                --rdrp-merged {input.merged} \
                --output      {output.itol} \
                2> {log.err}
            """

    rule itol_connection_by_rank:
        input:
            lookup            = _OUTDIR + "/6_fasttree/filtered/my_tree.renamed_lookup.tsv",
            esvirtu_table     = ancient(_ESVIRTU_INFO_TSV),
            known_known_dir   = ancient(_SPEARMAN_DIR + "/known_known_pair"),
            known_unknown_dir = ancient(_SPEARMAN_DIR + "/known_unknown_pair"),
            script            = "scripts/07_phylogenetic_tree/itol_connection.py",
        params:
            min_seed_nonzero   = _ITOL_MIN_SEED_NONZERO,
            min_target_nonzero = _ITOL_MIN_TARGET_NONZERO,
        output:
            expand(_OUTDIR + "/6_fasttree/filtered/itol/connection_by_rank/itol_connection_r{t}.txt",
                   t=["0.6", "0.7", "0.8", "0.9"]),
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered/itol_connection_by_rank.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output[0]})
            python {input.script} \
                --esvirtu-table        {input.esvirtu_table} \
                --lookup               {input.lookup} \
                --known-known-dir      {input.known_known_dir} \
                --known-unknown-dir    {input.known_unknown_dir} \
                --min-seed-nonzero     {params.min_seed_nonzero} \
                --min-target-nonzero   {params.min_target_nonzero} \
                --outdir               $(dirname {output[0]}) \
                2> {log.err}
            """

    rule itol_connection_by_virus:
        input:
            lookup            = _OUTDIR + "/6_fasttree/filtered/my_tree.renamed_lookup.tsv",
            esvirtu_table     = ancient(_ESVIRTU_INFO_TSV),
            known_known_dir   = ancient(_SPEARMAN_DIR + "/known_known_pair"),
            known_unknown_dir = ancient(_SPEARMAN_DIR + "/known_unknown_pair"),
            script            = "scripts/07_phylogenetic_tree/itol_connection_by_virus.py",
        params:
            min_seed_nonzero   = _ITOL_MIN_SEED_NONZERO,
            min_target_nonzero = _ITOL_MIN_TARGET_NONZERO,
        output:
            expand(_OUTDIR + "/6_fasttree/filtered/itol/connection_by_virus/itol_connection_r{t}.txt",
                   t=["0.6", "0.7", "0.8", "0.9"]),
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered/itol_connection_by_virus.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output[0]})
            python {input.script} \
                --esvirtu-table        {input.esvirtu_table} \
                --lookup               {input.lookup} \
                --known-known-dir      {input.known_known_dir} \
                --known-unknown-dir    {input.known_unknown_dir} \
                --min-seed-nonzero     {params.min_seed_nonzero} \
                --min-target-nonzero   {params.min_target_nonzero} \
                --outdir               $(dirname {output[0]}) \
                2> {log.err}
            """

    rule itol_gappa_phylum_colorstrip:
        input:
            phylum_tsv = _OUTDIR + "/5_gappa/filtered/phylum_annotation.tsv",
            script     = "scripts/07_phylogenetic_tree/itol_gappa_phylum_colorstrip.py",
        output:
            txt = _OUTDIR + "/6_fasttree/filtered/itol/colorstrip/itol_phylum_gappa.txt",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered/itol_phylum_gappa.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output.txt}) $(dirname {log.err})
            python {input.script} \
                --phylum-tsv {input.phylum_tsv} \
                --out        {output.txt} \
                2> {log.err}
            """


# ── FastTree filtered + ICTV ─────────────────────────────────────────────────
if _INCLUDE_FASTTREE_FILTERED_ICTV:
    rule concat_query_rt_filtered_ictv:
        input:
            query = _OUTDIR + "/3_split/filtered_ictv/query.fasta",
            rt    = _OUTDIR + "/6_fasttree/rt_outgroup.faa",
        output:
            faa = _OUTDIR + "/6_fasttree/filtered_ictv/query_with_rt.faa",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered_ictv/concat_query_rt.err",
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output.faa}) $(dirname {log.err})
            cat {input.query} {input.rt} > {output.faa} 2> {log.err}
            echo "Total sequences (filtered+ICTV): $(grep -c '^>' {output.faa})" >> {log.err}
            """

    rule fasttree_build_filtered_ictv:
        input:
            faa = _OUTDIR + "/6_fasttree/filtered_ictv/query_with_rt.faa",
        output:
            nwk = _OUTDIR + "/6_fasttree/filtered_ictv/my_tree.nwk",
        log:
            log = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered_ictv/fasttree.log",
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered_ictv/fasttree.err",
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
            FastTree -wag -gamma {_FASTTREE_BOOT_FLAG} \
                {input.faa} \
                > {output.nwk} \
                2> {log.log}
            """

    rule rename_tree_tips_filtered_ictv:
        input:
            nwk = _OUTDIR + "/6_fasttree/filtered_ictv/my_tree.nwk",
        output:
            nwk    = _OUTDIR + "/6_fasttree/filtered_ictv/my_tree.renamed.nwk",
            lookup = _OUTDIR + "/6_fasttree/filtered_ictv/my_tree.renamed_lookup.tsv",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered_ictv/rename.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {log.err})
            python scripts/07_phylogenetic_tree/rename_tree_tips.py \
                --input  {input.nwk} \
                --output {output.nwk} \
                --lookup {output.lookup} \
                2> {log.err}
            """

    rule classify_contigs_by_tree_filtered_ictv:
        input:
            tree   = _OUTDIR + "/6_fasttree/filtered_ictv/my_tree.renamed.nwk",
            lookup = _OUTDIR + "/6_fasttree/filtered_ictv/my_tree.renamed_lookup.tsv",
            vmr    = ancient(_ICTV_VMR_XLSX),
            script = "scripts/07_phylogenetic_tree/classify_contigs_by_tree.py",
        output:
            tsv = _OUTDIR + "/6_fasttree/filtered_ictv/contig_classification/contig_classification.tsv",
        params:
            vmr_sheet = _ICTV_VMR_SHEET,
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered_ictv/classify_contigs.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output.tsv}) $(dirname {log.err})
            python {input.script} \
                --tree      {input.tree} \
                --vmr       {input.vmr} \
                --vmr-sheet "{params.vmr_sheet}" \
                --lookup    {input.lookup} \
                --outdir    $(dirname {output.tsv}) \
                2> {log.err}
            """

    rule itol_ictv_colorstrip_filtered_ictv:
        input:
            lookup = _OUTDIR + "/6_fasttree/filtered_ictv/my_tree.renamed_lookup.tsv",
            vmr    = ancient(_ICTV_VMR_XLSX),
            script = "scripts/07_phylogenetic_tree/itol_ictv_colorstrip.py",
        output:
            expand(_OUTDIR + "/6_fasttree/filtered_ictv/itol/colorstrip/ictv/itol_ictv_{rank}.txt",
                   rank=["phylum", "class", "order", "family", "genus"]),
        params:
            vmr_sheet = _ICTV_VMR_SHEET,
            outdir    = _OUTDIR + "/6_fasttree/filtered_ictv/itol/colorstrip/ictv",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/filtered_ictv/itol_ictv_colorstrip.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p {params.outdir} $(dirname {log.err})
            python {input.script} \
                --vmr       {input.vmr} \
                --vmr-sheet "{params.vmr_sheet}" \
                --lookup    {input.lookup} \
                --outdir    {params.outdir} \
                2> {log.err}
            """


# ── FastTree unfiltered (all 4010+142 proteins, no contig-table filter) ──────

if _INCLUDE_FASTTREE_UNFILTERED:
    rule concat_query_rt_unfiltered:
        input:
            query = _OUTDIR + "/3_split/unfiltered/query.fasta",
            rt    = _OUTDIR + "/6_fasttree/rt_outgroup.faa",
        output:
            faa = _OUTDIR + "/6_fasttree/unfiltered/query_with_rt.faa",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/unfiltered/concat_query_rt.err",
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output.faa}) $(dirname {log.err})
            cat {input.query} {input.rt} > {output.faa} 2> {log.err}
            echo "Total sequences (full): $(grep -c '^>' {output.faa})" >> {log.err}
            """

    rule fasttree_build_unfiltered:
        input:
            faa = _OUTDIR + "/6_fasttree/unfiltered/query_with_rt.faa",
        output:
            nwk = _OUTDIR + "/6_fasttree/unfiltered/my_tree.nwk",
        log:
            log = "log/101_cluster_RDRP_tree_construction/6_fasttree/unfiltered/fasttree.log",
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/unfiltered/fasttree.err",
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
            FastTree -wag -gamma {_FASTTREE_BOOT_FLAG} \
                {input.faa} \
                > {output.nwk} \
                2> {log.log}
            """

    rule rename_tree_tips_unfiltered:
        input:
            nwk = _OUTDIR + "/6_fasttree/unfiltered/my_tree.nwk",
        output:
            nwk    = _OUTDIR + "/6_fasttree/unfiltered/my_tree.renamed.nwk",
            lookup = _OUTDIR + "/6_fasttree/unfiltered/my_tree.renamed_lookup.tsv",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/unfiltered/rename.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {log.err})
            python scripts/07_phylogenetic_tree/rename_tree_tips.py \
                --input  {input.nwk} \
                --output {output.nwk} \
                --lookup {output.lookup} \
                2> {log.err}
            """


# ── FastTree unfiltered + ICTV ────────────────────────────────────────────────

if _INCLUDE_FASTTREE_UNFILTERED_ICTV:
    rule concat_query_rt_unfiltered_ictv:
        input:
            query = _OUTDIR + "/3_split/unfiltered_ictv/query.fasta",
            rt    = _OUTDIR + "/6_fasttree/rt_outgroup.faa",
        output:
            faa = _OUTDIR + "/6_fasttree/unfiltered_ictv/query_with_rt.faa",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/unfiltered_ictv/concat_query_rt.err",
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {output.faa}) $(dirname {log.err})
            cat {input.query} {input.rt} > {output.faa} 2> {log.err}
            echo "Total sequences (unfiltered+ICTV): $(grep -c '^>' {output.faa})" >> {log.err}
            """

    rule fasttree_build_unfiltered_ictv:
        input:
            faa = _OUTDIR + "/6_fasttree/unfiltered_ictv/query_with_rt.faa",
        output:
            nwk = _OUTDIR + "/6_fasttree/unfiltered_ictv/my_tree.nwk",
        log:
            log = "log/101_cluster_RDRP_tree_construction/6_fasttree/unfiltered_ictv/fasttree.log",
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/unfiltered_ictv/fasttree.err",
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
            FastTree -wag -gamma {_FASTTREE_BOOT_FLAG} \
                {input.faa} \
                > {output.nwk} \
                2> {log.log}
            """

    rule rename_tree_tips_unfiltered_ictv:
        input:
            nwk = _OUTDIR + "/6_fasttree/unfiltered_ictv/my_tree.nwk",
        output:
            nwk    = _OUTDIR + "/6_fasttree/unfiltered_ictv/my_tree.renamed.nwk",
            lookup = _OUTDIR + "/6_fasttree/unfiltered_ictv/my_tree.renamed_lookup.tsv",
        log:
            err = "log/101_cluster_RDRP_tree_construction/6_fasttree/unfiltered_ictv/rename.err",
        conda:
            "../envs/python.yaml"
        threads: _SJ["threads"]
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = _SJ["runtime"],
            cpus_per_task   = _SJ["threads"],
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            """
            mkdir -p $(dirname {log.err})
            python scripts/07_phylogenetic_tree/rename_tree_tips.py \
                --input  {input.nwk} \
                --output {output.nwk} \
                --lookup {output.lookup} \
                2> {log.err}
            """


# ── Copy taxonomy annotations to 99_tables ────────────────────────────────────

if _EPANG_BATCHES:
    rule cp_gappa_phylum_annotation:
        wildcard_constraints:
            batch = "|".join(_EPANG_BATCHES.keys()),
        input:
            _OUTDIR + "/5_gappa/{batch}/phylum_annotation.tsv",
        output:
            _CP_TABLE_DIR + "/5_gappa/{batch}/phylum_annotation.tsv",
        threads: 1
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = 30,
            cpus_per_task   = 1,
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            "mkdir -p $(dirname {output}) && cp {input} {output}"

    rule cp_gappa_full_taxonomy_annotation:
        wildcard_constraints:
            batch = "|".join(_EPANG_BATCHES.keys()),
        input:
            _OUTDIR + "/5_gappa/{batch}/full_taxonomy_annotation.tsv",
        output:
            _CP_TABLE_DIR + "/5_gappa/{batch}/full_taxonomy_annotation.tsv",
        threads: 1
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = 30,
            cpus_per_task   = 1,
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            "mkdir -p $(dirname {output}) && cp {input} {output}"

if _INCLUDE_FASTTREE_FILTERED_ICTV:
    rule cp_contig_classification:
        input:
            _OUTDIR + "/6_fasttree/filtered_ictv/contig_classification/contig_classification.tsv",
        output:
            _CP_TABLE_DIR + "/6_fasttree/filtered_ictv/contig_classification.tsv",
        threads: 1
        resources:
            mem_mb_per_cpu  = _SJ["memory"],
            runtime         = 30,
            cpus_per_task   = 1,
            slurm_partition = _SJ["partition"],
            slurm_account   = _SJ["account"],
        shell:
            "mkdir -p $(dirname {output}) && cp {input} {output}"
