# RdRP phylogenetic PLACEMENT pipeline (EPA-ng)
#
# Places candidate RdRp query sequences onto the FIXED Neri et al. 2022 RCR90
# reference tree — this is placement, not de-novo tree search. The reference
# tree/MSA (config["epa_ng"]["ref_tree"] / ["ref_msa"]) must never be modified
# or re-optimized at any stage.
#
# Two query batches, same fixed reference:
#   batch "tier1ab"  = tier1a (ABCD canonical) + tier1b (CABD depermuted)
#   batch "tier1ab2" = tier1a + tier1b + tier2 (3-motif partial)
#
# Pipeline:
#   1. hmmbuild profile HMM from the reference RCR90 MSA (reference columns
#      preserved exactly — this defines the placement coordinate system)
#   2. hmmalign --trim queries into reference column coordinates (queries are
#      NOT added as new columns; output column count == reference MSA column count)
#   3. split hmmalign output back into reference-MSA-shaped and query-MSA-shaped
#      FASTA (EPA-ng requires the two as separate files, both already in the
#      same aligned column space)
#   4. epa-ng: place queries onto the fixed reference tree -> .jplace
#   5. gappa: graft placements onto the tree for visualization, and assign
#      taxonomy/lineage per query for Tier-1 ICTV-based reporting
#
# IMPORTANT — before running at scale on DelftBlue:
#   Verify every flag below against the ACTUALLY INSTALLED versions:
#     conda activate <env-from-envs/epang_gappa.yaml>  (or snakemake --use-conda)
#     hmmbuild -h
#     hmmalign -h
#     epa-ng --help
#     gappa examine graft --help
#     gappa examine assign --help
#   epa-ng's hmmer-alignment support and exact split-input flags have changed
#   across versions — do not trust flags here without checking --help on your
#   installed build. Flags marked PLACEHOLDER below are best-guess from the
#   documented stable CLI and must be confirmed before a full-scale run.

import os

configfile: "config/config.yaml"

_MOTIF_RESULTS      = "result/03_RDRP_identification/4_motif_search/motif_results"
_TIER1A_FAA         = _MOTIF_RESULTS + "/tier1a_ABCD.faa"
_TIER1B_FAA         = _MOTIF_RESULTS + "/tier1b_CABD_depermuted.faa"
_TIER2_FAA          = _MOTIF_RESULTS + "/tier2_3motif.faa"

# ICTV reference candidates go through the same motif pipeline (workflow 03,
# ICTV rules) and produce the same tier1a/1b/2 outputs — include them in the
# query batches so known ICTV taxa are placed onto the tree alongside contigs,
# giving a direct visual/taxonomic anchor for the Yangtze candidates.
_ICTV_MOTIF_RESULTS = "result/03_RDRP_identification/4_motif_search/ICTV/motif_results"
_ICTV_TIER1A_FAA    = _ICTV_MOTIF_RESULTS + "/tier1a_ABCD.faa"
_ICTV_TIER1B_FAA    = _ICTV_MOTIF_RESULTS + "/tier1b_CABD_depermuted.faa"
_ICTV_TIER2_FAA     = _ICTV_MOTIF_RESULTS + "/tier2_3motif.faa"

# ESvirtu reference candidates (workflow 04_03) also go through the same motif
# pipeline and produce the same tier1a/1b/2 outputs.
_ESV_MOTIF_RESULTS  = "result/04_03_esvirtu_rdrp_identification/7_motif_search/motif_results"
_ESV_TIER1A_FAA     = _ESV_MOTIF_RESULTS + "/tier1a_ABCD.faa"
_ESV_TIER1B_FAA     = _ESV_MOTIF_RESULTS + "/tier1b_CABD_depermuted.faa"
_ESV_TIER2_FAA      = _ESV_MOTIF_RESULTS + "/tier2_3motif.faa"

_OUTDIR  = "result/07_RDRP_tree_construction"
_REF_MSA      = config["epa_ng"]["ref_msa"]
_REF_TREE     = config["epa_ng"]["ref_tree"]
_MODEL        = config["epa_ng"]["model"]
_RVMT_INFO_TSV = config["epa_ng"]["rvmt_info_tsv"]
_TAXON_FILE   = _OUTDIR + "/0_taxon/RCR90_taxonomy.tsv"

# Whether to fold ICTV / ESvirtu reference candidates into the query batches
# alongside contig-derived candidates. Both are independently selectable via
# config["epa_ng"]["include_ictv"] and config["epa_ng"]["include_esvirtu"].
_INCLUDE_ICTV    = config["epa_ng"].get("include_ictv", False)
_INCLUDE_ESVIRTU = config["epa_ng"].get("include_esvirtu", False)
_INCLUDE_FASTTREE = config["epa_ng"].get("include_fasttree", False)

# Tier FAA files consumed by this workflow, staged into 0_inputdata/ as a frozen,
# self-contained snapshot (Step 0 / rule stage_inputs). Staged under
# "contig_"/"ictv_"/"esvirtu_" prefixes since the three sources share the same
# basenames (tier1a_ABCD.faa etc.) in their own output directories.
# Only the sources actually selected are staged/consumed.
_STAGED_INPUTS = {
    "contig_tier1a_ABCD.faa":            _TIER1A_FAA,
    "contig_tier1b_CABD_depermuted.faa": _TIER1B_FAA,
    "contig_tier2_3motif.faa":           _TIER2_FAA,
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

# batch_name -> list of STAGED input FAA files to concatenate as the query set
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
        # Step 0 — taxonomy file + staged FAA copies
        _TAXON_FILE,
        expand(_OUTDIR + "/0_inputdata/{f}", f=list(_STAGED_INPUTS.keys())),
        # Step 1 — reference HMM
        _OUTDIR + "/1_ref_hmm/RCR90.hmm",
        # Step 2-3 — per-batch aligned + split query/ref MSA
        expand(_OUTDIR + "/3_split/{batch}/query.fasta",     batch=list(_BATCHES.keys())),
        expand(_OUTDIR + "/3_split/{batch}/reference.fasta", batch=list(_BATCHES.keys())),
        # Step 4 — placements
        expand(_OUTDIR + "/4_epang/{batch}/epa_result.jplace", batch=list(_BATCHES.keys())),
        # Step 5a — gappa graft
        expand(_OUTDIR + "/5_gappa/{batch}/graft.newick", batch=list(_BATCHES.keys())),
        # Step 5b — gappa assign (taxonomy per query)
        expand(_OUTDIR + "/5_gappa/{batch}/per_query.tsv", batch=list(_BATCHES.keys())),
        # FastTree de-novo tree (contig + ESvirtu queries + RT outgroup)
        # Only requested when include_fasttree: true in config.
        *([
            _OUTDIR + "/6_fasttree/rt_outgroup.faa",
            _OUTDIR + "/6_fasttree/query_with_rt.faa",
            _OUTDIR + "/6_fasttree/my_tree.nwk",
            _OUTDIR + "/6_fasttree/itol_source.txt",
            _OUTDIR + "/6_fasttree/itol_tier.txt",
        ] if _INCLUDE_FASTTREE else []),


# ── Step 0: stage workflow 03 FAA files used by this workflow into 0_inputdata ─
# Copies (not symlinks) so 0_inputdata is a frozen, self-contained record of
# exactly which tier1a/1b/2 (contig + ICTV) sequences fed this placement run,
# independent of whatever workflow 03 produces afterward.

rule stage_inputs:
    input:
        f = lambda wc: _STAGED_INPUTS[wc.f],
    output:
        _OUTDIR + "/0_inputdata/{f}",
    log:
        err = "log/07_RDRP_tree_construction/0_inputdata/{f}.err",
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
# Reference alignment columns must be preserved exactly: hmmbuild does not
# alter the input MSA, it only builds a search/alignment profile from it.
# By default hmmbuild marks a column as a match state only if >50% of
# sequences have a residue there (--symfrac 0.5), so with a large, gappy
# reference MSA most columns can be dropped to insert state and the HMM ends
# up with far fewer match states than the MSA has columns (e.g. 416 of 874).
# That breaks the fixed reference-column coordinate system EPA-ng/gappa graft
# require. --symfrac 0.0 forces every alignment column to be a match state,
# so the HMM's match-state count always equals the reference MSA's column
# count exactly.

rule hmmbuild_reference:
    input:
        ref_msa = _REF_MSA,
    output:
        hmm = _OUTDIR + "/1_ref_hmm/RCR90.hmm",
    log:
        log = "log/07_RDRP_tree_construction/1_ref_hmm/hmmbuild.log",
        err = "log/07_RDRP_tree_construction/1_ref_hmm/hmmbuild.err",
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
# --trim removes insert states relative to the reference (match-state) columns,
# so the output stays in the reference's coordinate system instead of growing
# new columns for query-only insertions. CONFIRM: on your installed hmmer this
# may emit Stockholm by default — check `hmmalign -h` for the actual output
# format flag (commonly --outformat) before trusting the .sto path below.

rule hmmalign_query:
    input:
        hmm   = _OUTDIR + "/1_ref_hmm/RCR90.hmm",
        query = lambda wc: _BATCHES[wc.batch],
    output:
        sto = _OUTDIR + "/2_hmmalign/{batch}/aligned.sto",
    params:
        merged_query = _OUTDIR + "/2_hmmalign/{batch}/query_merged.faa",
    log:
        log = "log/07_RDRP_tree_construction/2_hmmalign/{batch}.log",
        err = "log/07_RDRP_tree_construction/2_hmmalign/{batch}.err",
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
# EPA-ng wants the reference MSA (unchanged, original columns) and query MSA
# (same columns, post --trim) as two separate FASTA files. We split by sequence
# ID: anything whose ID matches a reference-MSA sequence ID goes to
# reference.fasta, everything else (the queries we just aligned) goes to
# query.fasta. Both retain the identical (trimmed) column coordinate system,
# which is what EPA-ng requires for placement.
#
# hmmalign --trim only strips residues outside the first/last match state — it
# does NOT strip insert-state columns between match states, so the raw query
# Stockholm has more columns than the reference MSA. esl-alimask --rf-is-mask
# uses the "#=GC RF" annotation (x = match, . = insert) to keep only
# match-state columns, bringing query.fasta back to the same column count as
# reference.fasta (the reference MSA itself is already match-state-only, so it
# just needs reformatting to FASTA, no masking).

rule split_ref_query:
    input:
        sto     = _OUTDIR + "/2_hmmalign/{batch}/aligned.sto",
        ref_msa = _REF_MSA,
    output:
        ref_fasta   = _OUTDIR + "/3_split/{batch}/reference.fasta",
        query_fasta = _OUTDIR + "/3_split/{batch}/query.fasta",
    log:
        err = "log/07_RDRP_tree_construction/3_split/{batch}.err",
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
# CONFIRM all flags against `epa-ng --help` on the installed version before a
# full-scale run. The flags below reflect the documented stable EPA-ng CLI:
#   --tree            reference tree (Newick) — fixed, never re-optimized
#   --ref-msa         reference MSA (original columns)
#   --query           query MSA (same columns as ref-msa, post hmmalign --trim)
#   --model             RAxML-style model string (e.g. "WAG+G"), confirmed flag
#                       name from `epa-ng --help` on the installed version
#   --outdir          output directory (epa-ng writes epa_result.jplace there)
#   -T / --threads    thread count

rule epa_ng_place:
    input:
        tree     = _REF_TREE,
        ref_msa  = _OUTDIR + "/3_split/{batch}/reference.fasta",
        query    = _OUTDIR + "/3_split/{batch}/query.fasta",
    output:
        jplace = _OUTDIR + "/4_epang/{batch}/epa_result.jplace",
    params:
        outdir = _OUTDIR + "/4_epang/{batch}",
        model  = _MODEL,
    log:
        log = "log/07_RDRP_tree_construction/4_epang/{batch}.log",
        err = "log/07_RDRP_tree_construction/4_epang/{batch}.err",
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


# ── Step 5a: gappa graft — visualize placements on the reference tree ────────

rule gappa_graft:
    input:
        jplace = _OUTDIR + "/4_epang/{batch}/epa_result.jplace",
    output:
        newick = _OUTDIR + "/5_gappa/{batch}/graft.newick",
    params:
        outdir = _OUTDIR + "/5_gappa/{batch}",
    log:
        log = "log/07_RDRP_tree_construction/5_gappa/{batch}.graft.log",
        err = "log/07_RDRP_tree_construction/5_gappa/{batch}.graft.err",
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

# ── Prepare taxon file: RCR90 tip label -> semicolon-separated lineage ────────
# RiboV1.6_Info.tsv columns used:
#   col3  (RCR90)  = tree tip label, e.g. Rv4_170939
#   col12 (AfLvl)  = assignment level — only "Lvl 0 - Megatree leaves" are
#                    actual tree tips; Lvl 1+ are BLASTp-assigned and not in
#                    the reference tree, so they must be excluded
#   col16 (Phylum) col17 (Class) col18 (Order) col19 (Family) col20 (Genus)
# Output: two-column TSV (no header): tip_label \t Phylum;Class;Order;Family;Genus
# Rows where RCR90 is empty, starts with "rt." (outgroup), or AfLvl != Lvl 0 are skipped.

rule prepare_taxon_file:
    input:
        tsv = _RVMT_INFO_TSV,
    output:
        tsv = _TAXON_FILE,
    log:
        err = "log/07_RDRP_tree_construction/0_taxon/prepare_taxon_file.err",
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


# ── Step 5b: gappa assign — per-query taxonomy from placement ─────────────────

rule gappa_assign:
    input:
        jplace     = _OUTDIR + "/4_epang/{batch}/epa_result.jplace",
        taxon_file = _TAXON_FILE,
    output:
        tsv = _OUTDIR + "/5_gappa/{batch}/per_query.tsv",
    params:
        outdir = _OUTDIR + "/5_gappa/{batch}",
    log:
        log = "log/07_RDRP_tree_construction/5_gappa/{batch}.assign.log",
        err = "log/07_RDRP_tree_construction/5_gappa/{batch}.assign.err",
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
    # ════════════════════════════════════════════════════════════════════════════
    # Steps 6a-6c: de-novo FastTree (contig + ESvirtu queries + RT outgroup)
    #
    # Activated only when config["epa_ng"]["include_fasttree"] = true.
    #
    # Methodology follows Neri et al. 2022:
    #   - WAG substitution model, gamma rate variation (FastTree -wag -gamma)
    #   - All sequences already in the same 874-column coordinate system
    #     (produced by hmmalign --trim + esl-alimask via Steps 2-3 above)
    #   - 10 reverse-transcriptase outgroup sequences extracted directly from
    #     the reference MSA: RT sequences have IDs starting with "rt." (e.g.
    #     rt.489297600). Extracted in their already-aligned 874-column form —
    #     no re-alignment needed.
    #   - Query batch used: tier1ab2 (tier1a + tier1b + tier2) — change the
    #     hard-coded batch name below if you want a different subset.
    #
    # Output: result/07_RDRP_tree_construction/6_fasttree/my_tree.nwk
    # ════════════════════════════════════════════════════════════════════════════

    # ── Step 6a: extract RT outgroup sequences from reference MSA ───────────────
    # RT sequences have IDs starting with "rt." (e.g. rt.489297600).
    # Extracted in their already-aligned form (esl-reformat afa converts the
    # reference MSA to FASTA, then seqkit grep pulls the rt.* headers).

    rule extract_rt_outgroup:
        input:
            ref_msa = _REF_MSA,
        output:
            faa = _OUTDIR + "/6_fasttree/rt_outgroup.faa",
        log:
            err = "log/07_RDRP_tree_construction/6_fasttree/extract_rt_outgroup.err",
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
            # Convert reference MSA to FASTA, then keep only rt.* sequences (RTs)
            esl-reformat afa {input.ref_msa} 2>> {log.err} \
                | seqkit grep --use-regexp --pattern '^rt\\.' --threads {threads} \
                > {output.faa} 2>> {log.err}
            echo "RT outgroup count: $(grep -c '^>' {output.faa})" >> {log.err}
            """


    # ── Step 6b: concatenate query + RT outgroup ─────────────────────────────────
    # Uses the tier1ab2 query batch (already split in Step 3). All sequences are
    # in the same 874-column space — no re-alignment needed.

    rule concat_query_rt:
        input:
            query = _OUTDIR + "/3_split/tier1ab2/query.fasta",
            rt    = _OUTDIR + "/6_fasttree/rt_outgroup.faa",
        output:
            faa = _OUTDIR + "/6_fasttree/query_with_rt.faa",
        log:
            err = "log/07_RDRP_tree_construction/6_fasttree/concat_query_rt.err",
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


    # ── Step 6c: FastTree de-novo tree (WAG + gamma, Neri 2022 parameters) ──────

    rule fasttree_build:
        input:
            faa = _OUTDIR + "/6_fasttree/query_with_rt.faa",
        output:
            nwk = _OUTDIR + "/6_fasttree/my_tree.nwk",
        log:
            log = "log/07_RDRP_tree_construction/6_fasttree/fasttree.log",
            err = "log/07_RDRP_tree_construction/6_fasttree/fasttree.err",
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


    # ── Step 6d: iTOL color-strip decoration files ───────────────────────────────
    # Produces two files:
    #   itol_source.txt — sequence source (rt / esvirtu / assembled)
    #   itol_tier.txt   — RdRp tier (tier1a / tier1b / tier2 / n/a)

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
            err = "log/07_RDRP_tree_construction/6_fasttree/itol_colorstrip.err",
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
