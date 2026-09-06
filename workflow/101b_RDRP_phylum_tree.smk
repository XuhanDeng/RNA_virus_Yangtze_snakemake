# Per-Phylum/Class/Order/Family RdRp tree construction.
# Query proteins ({taxon}.faa) are merged per mode/rank/taxon from THREE
# pipelines (mode in full_length/palm_core/palm_extended):
#   03  result/03_RDRP_identification/22_phylum_cluster/{mode}_{rank}
#   03b result/03b_esvirtu_rdrp_identification/13_phylum_cluster/{mode}_{rank}
#   03c result/03c_ref_palm_annot/2_rdrp_region/{rank}/{taxon}/{03c_file}
#       (03c_file: rdrp_full.faa / palm_core.faa / palm_extended.faa)
# Reference/outgroup sequences all come from 03c's RVMT_ref subfolder (03's
# own {taxon}_ref.faa is NOT staged separately -- it's the same underlying
# sequences 03c already re-emits via rdrp_full.faa, just before palm_annot
# added motif tags -- see rule stage_phylum_inputs).
# Staging (Step 0) runs unconditionally for all 3 modes, independent of which
# tree-building mode(s) are toggled on below -- see rule stage_phylum_inputs.
# Run after 03_RDRP_identification.smk (Steps 20-22), 03b, and 03c are complete.

import os
import types

configfile: "config/config.yaml"

_RANKS   = ["phylum", "class", "order", "family"]
_MODES   = ["full_length", "palm_core", "palm_extended"]
_SRC_DIR      = "result/03_RDRP_identification/22_phylum_cluster"
_SRC_DIR_03B  = "result/03b_esvirtu_rdrp_identification/13_phylum_cluster"
_SRC_DIR_03C  = "result/03c_ref_palm_annot/2_rdrp_region"
_OUTDIR  = "result/101b_RDRP_phylum_tree"

# 03c's per-taxon region file name, keyed by mode (same naming used in 03c's
# extract_ref_palm_regions rule output).
_MODE_03C_FILE = {
    "full_length":   "rdrp_full.faa",
    "palm_core":     "palm_core.faa",
    "palm_extended": "palm_extended.faa",
}


# rdrp_tree.{full_length,palm_core,palm_extended} each toggle their own
# tree-building block further below. full_length (this file's original
# implementation) aligns full-length RdRp proteins directly, which was found
# to blow up alignment width badly (8,800-14,000+ columns) for phyla with
# polyprotein/fusion-length outliers -- see conversation notes.
# palm_core/palm_extended avoid this by aligning a motif-bounded region
# instead. Any combination of the 3 modes may be enabled at once via config.
#
# rule all (the single Snakemake no-target default, defined at the very top
# of the file) collects whichever mode(s) are actually enabled -- built as a
# plain Python list below rather than three separate `rule all_<mode>`
# targets, so a bare `snakemake` run with no target builds every
# config-enabled mode in one pass.

_ALL_TARGETS = []

if config["rdrp_tree"]["full_length"]:
    _ALL_TARGETS += [
        *expand(_OUTDIR + "/0_input/full_length/{rank}", rank=_RANKS),
        *expand(_OUTDIR + "/0_input/full_length/{rank}_taxa.txt", rank=_RANKS),
    ]

if config["rdrp_tree"]["palm_core"]:
    _ALL_TARGETS += [
        *expand(_OUTDIR + "/0_input/palm_core/{rank}", rank=_RANKS),
        *expand(_OUTDIR + "/0_input/palm_core/{rank}_taxa.txt", rank=_RANKS),
    ]

if config["rdrp_tree"]["palm_extended"]:
    _ALL_TARGETS += [
        *expand(_OUTDIR + "/0_input/palm_extended/{rank}", rank=_RANKS),
        *expand(_OUTDIR + "/0_input/palm_extended/{rank}_taxa.txt", rank=_RANKS),
    ]


def all_enabled_trimmed_and_trees(wildcards):
    files = []
    if config["rdrp_tree"]["full_length"]:
        for rank in ["phylum", "class"]:
            files += all_alignments(types.SimpleNamespace(rank=rank))
            files += all_trees(types.SimpleNamespace(rank=rank))
    if config["rdrp_tree"]["palm_core"]:
        for rank in ["phylum", "class"]:
            files += all_alignments_palm_core(types.SimpleNamespace(rank=rank))
            files += all_trees_palm_core(types.SimpleNamespace(rank=rank))
    if config["rdrp_tree"]["palm_extended"]:
        for rank in ["phylum", "class"]:
            files += all_alignments_palm_extended(types.SimpleNamespace(rank=rank))
            files += all_trees_palm_extended(types.SimpleNamespace(rank=rank))
    return files


rule all:
    input:
        _ALL_TARGETS,
        all_enabled_trimmed_and_trees,


# ── Step 0: stage a self-contained snapshot of the per-taxon input fastas ────
# Runs for all 3 modes unconditionally (outside any rdrp_tree.* toggle) --
# staging is shared groundwork, not specific to which tree-building mode(s)
# are enabled below. No merging happens here: each source is copied as-is
# into its own subfolder --
#   0_input/{mode}/{rank}/new/{taxon}.faa       -- 03's assembled-contig proteins
#   0_input/{mode}/{rank}/esvirtu/{taxon}.faa   -- 03b's esvirtu-reference proteins
#   0_input/{mode}/{rank}/RVMT_ref/{taxon}.faa  -- 03c's region-extracted RVMT
#                                                   reference (renamed from
#                                                   {mode_03c_file}); this is
#                                                   the sole RVMT reference
#                                                   source for all 3 modes,
#                                                   including full_length --
#                                                   03's own {taxon}_ref.faa is
#                                                   NOT staged separately, to
#                                                   avoid double-counting the
#                                                   same sequences (see shell
#                                                   block below for detail)
# Merging the three query sources into one alignment input happens later, at
# combine_query_ref (per mode).

rule stage_phylum_inputs:
    input:
        srcdir_03  = _SRC_DIR + "/{mode}_{rank}",
        srcdir_03b = _SRC_DIR_03B + "/{mode}_{rank}",
    output:
        outdir = directory(_OUTDIR + "/0_input/{mode}/{rank}"),
    wildcard_constraints:
        mode = "full_length|palm_core|palm_extended",
        rank = "phylum|class|order|family",
    params:
        rdrpregion_dir = _SRC_DIR_03C + "/{rank}",
        mode_03c_file  = lambda wc: _MODE_03C_FILE[wc.mode],
    log:
        err = "log/101b_RDRP_phylum_tree/0_input/stage_phylum_inputs_{mode}_{rank}.err",
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        shopt -s nullglob
        mkdir -p {output.outdir}/new {output.outdir}/esvirtu {output.outdir}/RVMT_ref $(dirname {log.err})
        : > {log.err}

        # new: assembled-contig query proteins (03), as-is (excluding _ref.faa).
        for f in {input.srcdir_03}/*.faa; do
            [[ "$(basename "$f" .faa)" != *_ref ]] && cp "$f" {output.outdir}/new/
        done

        # esvirtu: esvirtu-reference query proteins (03b), as-is (excluding _ref.faa).
        for f in {input.srcdir_03b}/*.faa; do
            [[ "$(basename "$f" .faa)" != *_ref ]] && cp "$f" {output.outdir}/esvirtu/
        done

        # RVMT_ref: region-extracted RVMT reference (03c), renamed {{taxon}}.faa
        # (from {params.mode_03c_file}) per taxon subfolder. This is the sole
        # RVMT reference source for all 3 modes (including full_length) -- 03's
        # own {{taxon}}_ref.faa is NOT copied separately here, since it's the
        # same 879-ish sequences as RVMT_ref/{{taxon}}.faa (03c stages 03's
        # {{taxon}}_ref.faa as its own input, then re-emits it via rdrp_full.faa
        # with palm_annot motif tags added) -- copying both would double-count
        # every RVMT reference sequence in the alignment.
        for f in {params.rdrpregion_dir}/*/{params.mode_03c_file}; do
            taxon=$(basename "$(dirname "$f")")
            cp "$f" {output.outdir}/RVMT_ref/${{taxon}}.faa
        done
        true
        """


if config["rdrp_tree"]["full_length"]:

    # Reciprocal cross-Phylum outgroup: Duplornaviricota (dsRNA) is fairly
    # distinct from the other four (+ssRNA/-ssRNA) phyla, so each side is used
    # as the other's outgroup. N_OUTGROUP reference sequences are taken (in
    # file order, from that Phylum's {Phylum}_ref.faa) and added to the
    # Phylum-level alignment only -- Class/Order/Family sub-alignments are
    # subsets of an already-rooted Phylum tree and don't need their own outgroup.
    N_OUTGROUP = config["rdrp_tree"]["n_outgroup"]
    # The 5 canonical RNA virus phyla for the reciprocal Duplornaviricota-vs-rest
    # outgroup scheme. Only phyla actually present as ingroup taxa (i.e. that have
    # a {Phylum}.faa in 0_input/full_length/phylum/) get an outgroup file built --
    # this is resolved dynamically per-run via _outgroup_phylum_map() below, so a
    # future dataset that does/doesn't contain a given phylum needs no code change.
    _KNOWN_PHYLA = [
        "Duplornaviricota", "Pisuviricota", "Kitrinoviricota",
        "Lenarviricota", "Negarnaviricota",
    ]


    def _outgroup_phylum_map():
        """{ingroup phylum -> [candidate outgroup source phyla]}, restricted to
        phyla actually present in 0_input/full_length/phylum/RVMT_ref/."""
        indir = _OUTDIR + "/0_input/full_length/phylum/RVMT_ref"
        if not os.path.isdir(indir):
            return {}
        present = {
            p for p in _KNOWN_PHYLA
            if os.path.exists(f"{indir}/{p}.faa")
        }
        mapping = {}
        if "Duplornaviricota" in present:
            others = sorted(present - {"Duplornaviricota"})
            if others:
                mapping["Duplornaviricota"] = others
            for p in others:
                mapping[p] = ["Duplornaviricota"]
        return mapping


    # ── Step 1: enumerate taxa actually present per rank ──────────────────────
    # Taxon names (Pisuviricota, Marnaviridae, Unclassified, ...) are only
    # known after Step 0 (stage_phylum_inputs) has run. A Snakemake checkpoint
    # was used here previously, but (same as workflow 03/03c) it did not
    # reliably re-trigger downstream job discovery through rule all's
    # dependency chain on this Snakemake version -- see conversation notes.
    # list_taxa is a plain rule that writes the taxon names to a text file;
    # _taxa_for_rank reads that file directly at DAG-build time. This means
    # list_taxa's output must already exist on disk before the alignment/tree
    # targets show up in `rule all` -- rule all explicitly requests it first
    # below, so it always builds before anything else in this file.
    # Taxa are the union of whatever new/esvirtu/RVMT_ref staged for this
    # rank (each in its own subfolder under
    # 0_input/full_length/{rank}/{new,esvirtu,RVMT_ref}/).

    rule list_taxa:
        input:
            indir = _OUTDIR + "/0_input/full_length/{rank}",
        output:
            txt = _OUTDIR + "/0_input/full_length/{rank}_taxa.txt",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            shopt -s nullglob
            {{
                ls {input.indir}/new/*.faa {input.indir}/esvirtu/*.faa {input.indir}/RVMT_ref/*.faa 2> /dev/null \
                    | xargs -n1 basename 2> /dev/null | sed 's/\\.faa$//'
                true
            }} | sort -u > {output.txt}
            """


    def _taxa_for_rank(wc):
        taxa_file = _OUTDIR + f"/0_input/full_length/{wc.rank}_taxa.txt"
        if not os.path.exists(taxa_file):
            return []
        with open(taxa_file) as fh:
            return [line.strip() for line in fh if line.strip()]


    def combined_alignment_input(wc):
        indir = _OUTDIR + f"/0_input/full_length/{wc.rank}"
        files = []
        for src in ("new", "esvirtu", "RVMT_ref"):
            part = f"{indir}/{src}/{wc.taxon}.faa"
            if os.path.exists(part):
                files.append(part)
        if wc.rank == "phylum" and wc.taxon in _outgroup_phylum_map():
            files.append(_OUTDIR + f"/1_outgroup/{wc.taxon}.faa")
        return files


    def all_combined_fastas(wc):
        taxa = _taxa_for_rank(wc)
        return expand(
            _OUTDIR + "/1_combined/full_length_{rank}/{taxon}.faa",
            rank=wc.rank, taxon=taxa,
        )


    def all_alignments(wc):
        taxa = _taxa_for_rank(wc)
        return expand(
            _OUTDIR + "/3_trimmed/full_length_{rank}/{taxon}.trimmed.fasta",
            rank=wc.rank, taxon=taxa,
        )


    def all_trees(wc):
        taxa = _taxa_for_rank(wc)
        return expand(
            _OUTDIR + "/4_fasttree/full_length_{rank}/{taxon}.nwk",
            rank=wc.rank, taxon=taxa,
        )


    # ── Step 1b: extract N_OUTGROUP reference sequences per outgroup Phylum ──────
    # One outgroup file per ingroup Phylum (not per outgroup Phylum), since
    # Duplornaviricota draws its outgroup sequences from 4 different phyla, one
    # at a time, rather than pooling all of them together.
    # Source: 03c_ref_palm_annot's per-Phylum rdrp_full.faa (region-extracted
    # RVMT reference, full-length variant) rather than 03's raw {taxon}_ref.faa
    # -- same underlying reference sequences, already staged by 03c.

    def outgroup_ref_source(wc):
        outgroup_phyla = _outgroup_phylum_map()[wc.taxon]
        # Duplornaviricota case: multiple candidate outgroup phyla listed --
        # take the first one deterministically.
        return f"result/03c_ref_palm_annot/2_rdrp_region/phylum/{outgroup_phyla[0]}/rdrp_full.faa"

    rule extract_outgroup:
        input:
            ref = outgroup_ref_source,
        output:
            fasta = _OUTDIR + "/1_outgroup/{taxon}.faa",
        wildcard_constraints:
            taxon = "|".join(_KNOWN_PHYLA),
        conda:
            "../envs/mafft_hmmer_seqkit.yaml"
        log:
            err = "log/101b_RDRP_phylum_tree/1_outgroup/{taxon}.err",
        params:
            n = N_OUTGROUP,
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.fasta}) $(dirname {log.err})
            seqkit head -n {params.n} {input.ref} 2> {log.err} \
                | sed 's/^>/>OUTGROUP_/' > {output.fasta}
            """


    # ── Step 2: combine query + reference sequences per taxon ────────────────────

    rule combine_query_ref:
        input:
            files = combined_alignment_input,
        output:
            fasta = _OUTDIR + "/1_combined/full_length_{rank}/{taxon}.faa",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        log:
            err = "log/101b_RDRP_phylum_tree/1_combined/full_length_{rank}_{taxon}.err",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.fasta}) $(dirname {log.err})
            cat {input.files} > {output.fasta} 2> {log.err}
            """


    # ── Step 3: MAFFT alignment (--auto) ──────────────────────────────────────────

    rule mafft_align:
        input:
            fasta = _OUTDIR + "/1_combined/full_length_{rank}/{taxon}.faa",
        output:
            aln = _OUTDIR + "/2_aligned/full_length_{rank}/{taxon}.aln.fasta",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        conda:
            "../envs/mafft_hmmer_seqkit.yaml"
        log:
            err = "log/101b_RDRP_phylum_tree/2_aligned/full_length_{rank}_{taxon}.err",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.aln}) $(dirname {log.err})
            mafft --auto --thread {threads} {input.fasta} > {output.aln} 2> {log.err}
            """


    # ── Step 4: trimAl trimming ───────────────────────────────────────────────────

    rule trimal_trim:
        input:
            aln = _OUTDIR + "/2_aligned/full_length_{rank}/{taxon}.aln.fasta",
        output:
            trimmed = _OUTDIR + "/3_trimmed/full_length_{rank}/{taxon}.trimmed.fasta",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        conda:
            "../envs/trimal.yaml"
        log:
            err = "log/101b_RDRP_phylum_tree/3_trimmed/full_length_{rank}_{taxon}.err",
        params:
            resoverlap = config["rdrp_tree"]["trimal_resoverlap"],
            seqoverlap = config["rdrp_tree"]["trimal_seqoverlap"],
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.trimmed}) $(dirname {log.err})
            trimal -in {input.aln} -out {output.trimmed} \
                -resoverlap {params.resoverlap} -seqoverlap {params.seqoverlap} 2> {log.err}
            """


    # ── Step 5: FastTree maximum-likelihood tree (WAG + gamma) ───────────────────

    rule fasttree_ml:
        input:
            trimmed = _OUTDIR + "/3_trimmed/full_length_{rank}/{taxon}.trimmed.fasta",
        output:
            nwk = _OUTDIR + "/4_fasttree/full_length_{rank}/{taxon}.nwk",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        conda:
            "../envs/fasttree.yaml"
        log:
            log = "log/101b_RDRP_phylum_tree/4_fasttree/full_length_{rank}_{taxon}.log",
            err = "log/101b_RDRP_phylum_tree/4_fasttree/full_length_{rank}_{taxon}.err",
        params:
            flags = config["rdrp_tree"]["fasttree_flags"],
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.nwk}) $(dirname {log.log})
            export OMP_NUM_THREADS={threads}
            FastTree {params.flags} {input.trimmed} > {output.nwk} 2> {log.log}
            """


# rdrp_tree.palm_core toggles the block below. Same pipeline shape as the
# full_length block above, but aligning palm_core-region proteins (motif
# A/B/C-bounded cut) instead -- avoids the alignment-width blowup seen with
# full_length for phyla with polyprotein/fusion-length outliers. There is no
# {taxon}_ref.faa reference in 0_input/palm_core/{rank}/ (only full_length has
# one, from 03's Step 22) -- Phylum-level rooting relies on the outgroup file
# alone (03c's palm_core.faa, region-matched).
if config["rdrp_tree"]["palm_core"]:

    N_OUTGROUP = config["rdrp_tree"]["n_outgroup"]
    _KNOWN_PHYLA = [
        "Duplornaviricota", "Pisuviricota", "Kitrinoviricota",
        "Lenarviricota", "Negarnaviricota",
    ]


    def _outgroup_phylum_map_palm_core():
        """{ingroup phylum -> [candidate outgroup source phyla]}, restricted to
        phyla actually present (in any of new/esvirtu/RVMT_ref) in
        0_input/palm_core/phylum/."""
        indir = _OUTDIR + "/0_input/palm_core/phylum"
        if not os.path.isdir(indir):
            return {}
        present = {
            p for p in _KNOWN_PHYLA
            if any(
                os.path.exists(f"{indir}/{src}/{p}.faa")
                for src in ("new", "esvirtu", "RVMT_ref")
            )
        }
        mapping = {}
        if "Duplornaviricota" in present:
            others = sorted(present - {"Duplornaviricota"})
            if others:
                mapping["Duplornaviricota"] = others
            for p in others:
                mapping[p] = ["Duplornaviricota"]
        return mapping



    rule list_taxa_palm_core:
        input:
            indir = _OUTDIR + "/0_input/palm_core/{rank}",
        output:
            txt = _OUTDIR + "/0_input/palm_core/{rank}_taxa.txt",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            shopt -s nullglob
            {{
                ls {input.indir}/new/*.faa {input.indir}/esvirtu/*.faa {input.indir}/RVMT_ref/*.faa 2> /dev/null \
                    | xargs -n1 basename 2> /dev/null | sed 's/\\.faa$//'
                true
            }} | sort -u > {output.txt}
            """


    def _taxa_for_rank_palm_core(wc):
        taxa_file = _OUTDIR + f"/0_input/palm_core/{wc.rank}_taxa.txt"
        if not os.path.exists(taxa_file):
            return []
        with open(taxa_file) as fh:
            return [line.strip() for line in fh if line.strip()]


    def combined_alignment_input_palm_core(wc):
        indir = _OUTDIR + f"/0_input/palm_core/{wc.rank}"
        files = []
        for src in ("new", "esvirtu", "RVMT_ref"):
            part = f"{indir}/{src}/{wc.taxon}.faa"
            if os.path.exists(part):
                files.append(part)
        if wc.rank == "phylum" and wc.taxon in _outgroup_phylum_map_palm_core():
            files.append(_OUTDIR + f"/1_outgroup/palm_core/{wc.taxon}.faa")
        return files


    def all_alignments_palm_core(wc):
        taxa = _taxa_for_rank_palm_core(wc)
        return expand(
            _OUTDIR + "/3_trimmed/palm_core_{rank}/{taxon}.trimmed.fasta",
            rank=wc.rank, taxon=taxa,
        )


    def all_trees_palm_core(wc):
        taxa = _taxa_for_rank_palm_core(wc)
        return expand(
            _OUTDIR + "/4_fasttree/palm_core_{rank}/{taxon}.nwk",
            rank=wc.rank, taxon=taxa,
        )


    def outgroup_ref_source_palm_core(wc):
        outgroup_phyla = _outgroup_phylum_map_palm_core()[wc.taxon]
        return f"result/03c_ref_palm_annot/2_rdrp_region/phylum/{outgroup_phyla[0]}/palm_core.faa"

    rule extract_outgroup_palm_core:
        input:
            ref = outgroup_ref_source_palm_core,
        output:
            fasta = _OUTDIR + "/1_outgroup/palm_core/{taxon}.faa",
        wildcard_constraints:
            taxon = "|".join(_KNOWN_PHYLA),
        conda:
            "../envs/mafft_hmmer_seqkit.yaml"
        log:
            err = "log/101b_RDRP_phylum_tree/1_outgroup/palm_core_{taxon}.err",
        params:
            n = N_OUTGROUP,
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.fasta}) $(dirname {log.err})
            seqkit head -n {params.n} {input.ref} 2> {log.err} \
                | sed 's/^>/>OUTGROUP_/' > {output.fasta}
            """


    rule combine_query_ref_palm_core:
        input:
            files = combined_alignment_input_palm_core,
        output:
            fasta = _OUTDIR + "/1_combined/palm_core_{rank}/{taxon}.faa",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        log:
            err = "log/101b_RDRP_phylum_tree/1_combined/palm_core_{rank}_{taxon}.err",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.fasta}) $(dirname {log.err})
            cat {input.files} > {output.fasta} 2> {log.err}
            """


    rule mafft_align_palm_core:
        input:
            fasta = _OUTDIR + "/1_combined/palm_core_{rank}/{taxon}.faa",
        output:
            aln = _OUTDIR + "/2_aligned/palm_core_{rank}/{taxon}.aln.fasta",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        conda:
            "../envs/mafft_hmmer_seqkit.yaml"
        log:
            err = "log/101b_RDRP_phylum_tree/2_aligned/palm_core_{rank}_{taxon}.err",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.aln}) $(dirname {log.err})
            mafft --auto --thread {threads} {input.fasta} > {output.aln} 2> {log.err}
            """


    rule trimal_trim_palm_core:
        input:
            aln = _OUTDIR + "/2_aligned/palm_core_{rank}/{taxon}.aln.fasta",
        output:
            trimmed = _OUTDIR + "/3_trimmed/palm_core_{rank}/{taxon}.trimmed.fasta",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        conda:
            "../envs/trimal.yaml"
        log:
            err = "log/101b_RDRP_phylum_tree/3_trimmed/palm_core_{rank}_{taxon}.err",
        params:
            resoverlap = config["rdrp_tree"]["trimal_resoverlap"],
            seqoverlap = config["rdrp_tree"]["trimal_seqoverlap"],
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.trimmed}) $(dirname {log.err})
            trimal -in {input.aln} -out {output.trimmed} \
                -resoverlap {params.resoverlap} -seqoverlap {params.seqoverlap} 2> {log.err}
            """


    rule fasttree_ml_palm_core:
        input:
            trimmed = _OUTDIR + "/3_trimmed/palm_core_{rank}/{taxon}.trimmed.fasta",
        output:
            nwk = _OUTDIR + "/4_fasttree/palm_core_{rank}/{taxon}.nwk",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        conda:
            "../envs/fasttree.yaml"
        log:
            log = "log/101b_RDRP_phylum_tree/4_fasttree/palm_core_{rank}_{taxon}.log",
            err = "log/101b_RDRP_phylum_tree/4_fasttree/palm_core_{rank}_{taxon}.err",
        params:
            flags = config["rdrp_tree"]["fasttree_flags"],
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.nwk}) $(dirname {log.log})
            export OMP_NUM_THREADS={threads}
            FastTree {params.flags} {input.trimmed} > {output.nwk} 2> {log.log}
            """


# rdrp_tree.palm_extended toggles the block below. Same shape as palm_core
# above, but using the wider palm_extended region (core +/- flank aa).
if config["rdrp_tree"]["palm_extended"]:

    N_OUTGROUP = config["rdrp_tree"]["n_outgroup"]
    _KNOWN_PHYLA = [
        "Duplornaviricota", "Pisuviricota", "Kitrinoviricota",
        "Lenarviricota", "Negarnaviricota",
    ]


    def _outgroup_phylum_map_palm_extended():
        """{ingroup phylum -> [candidate outgroup source phyla]}, restricted to
        phyla actually present (in any of new/esvirtu/RVMT_ref) in
        0_input/palm_extended/phylum/."""
        indir = _OUTDIR + "/0_input/palm_extended/phylum"
        if not os.path.isdir(indir):
            return {}
        present = {
            p for p in _KNOWN_PHYLA
            if any(
                os.path.exists(f"{indir}/{src}/{p}.faa")
                for src in ("new", "esvirtu", "RVMT_ref")
            )
        }
        mapping = {}
        if "Duplornaviricota" in present:
            others = sorted(present - {"Duplornaviricota"})
            if others:
                mapping["Duplornaviricota"] = others
            for p in others:
                mapping[p] = ["Duplornaviricota"]
        return mapping


    rule list_taxa_palm_extended:
        input:
            indir = _OUTDIR + "/0_input/palm_extended/{rank}",
        output:
            txt = _OUTDIR + "/0_input/palm_extended/{rank}_taxa.txt",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            shopt -s nullglob
            {{
                ls {input.indir}/new/*.faa {input.indir}/esvirtu/*.faa {input.indir}/RVMT_ref/*.faa 2> /dev/null \
                    | xargs -n1 basename 2> /dev/null | sed 's/\\.faa$//'
                true
            }} | sort -u > {output.txt}
            """


    def _taxa_for_rank_palm_extended(wc):
        taxa_file = _OUTDIR + f"/0_input/palm_extended/{wc.rank}_taxa.txt"
        if not os.path.exists(taxa_file):
            return []
        with open(taxa_file) as fh:
            return [line.strip() for line in fh if line.strip()]


    def combined_alignment_input_palm_extended(wc):
        indir = _OUTDIR + f"/0_input/palm_extended/{wc.rank}"
        files = []
        for src in ("new", "esvirtu", "RVMT_ref"):
            part = f"{indir}/{src}/{wc.taxon}.faa"
            if os.path.exists(part):
                files.append(part)
        if wc.rank == "phylum" and wc.taxon in _outgroup_phylum_map_palm_extended():
            files.append(_OUTDIR + f"/1_outgroup/palm_extended/{wc.taxon}.faa")
        return files


    def all_alignments_palm_extended(wc):
        taxa = _taxa_for_rank_palm_extended(wc)
        return expand(
            _OUTDIR + "/3_trimmed/palm_extended_{rank}/{taxon}.trimmed.fasta",
            rank=wc.rank, taxon=taxa,
        )


    def all_trees_palm_extended(wc):
        taxa = _taxa_for_rank_palm_extended(wc)
        return expand(
            _OUTDIR + "/4_fasttree/palm_extended_{rank}/{taxon}.nwk",
            rank=wc.rank, taxon=taxa,
        )


    def outgroup_ref_source_palm_extended(wc):
        outgroup_phyla = _outgroup_phylum_map_palm_extended()[wc.taxon]
        return f"result/03c_ref_palm_annot/2_rdrp_region/phylum/{outgroup_phyla[0]}/palm_extended.faa"

    rule extract_outgroup_palm_extended:
        input:
            ref = outgroup_ref_source_palm_extended,
        output:
            fasta = _OUTDIR + "/1_outgroup/palm_extended/{taxon}.faa",
        wildcard_constraints:
            taxon = "|".join(_KNOWN_PHYLA),
        conda:
            "../envs/mafft_hmmer_seqkit.yaml"
        log:
            err = "log/101b_RDRP_phylum_tree/1_outgroup/palm_extended_{taxon}.err",
        params:
            n = N_OUTGROUP,
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.fasta}) $(dirname {log.err})
            seqkit head -n {params.n} {input.ref} 2> {log.err} \
                | sed 's/^>/>OUTGROUP_/' > {output.fasta}
            """


    rule combine_query_ref_palm_extended:
        input:
            files = combined_alignment_input_palm_extended,
        output:
            fasta = _OUTDIR + "/1_combined/palm_extended_{rank}/{taxon}.faa",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        log:
            err = "log/101b_RDRP_phylum_tree/1_combined/palm_extended_{rank}_{taxon}.err",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.fasta}) $(dirname {log.err})
            cat {input.files} > {output.fasta} 2> {log.err}
            """


    rule mafft_align_palm_extended:
        input:
            fasta = _OUTDIR + "/1_combined/palm_extended_{rank}/{taxon}.faa",
        output:
            aln = _OUTDIR + "/2_aligned/palm_extended_{rank}/{taxon}.aln.fasta",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        conda:
            "../envs/mafft_hmmer_seqkit.yaml"
        log:
            err = "log/101b_RDRP_phylum_tree/2_aligned/palm_extended_{rank}_{taxon}.err",
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.aln}) $(dirname {log.err})
            mafft --auto --thread {threads} {input.fasta} > {output.aln} 2> {log.err}
            """


    rule trimal_trim_palm_extended:
        input:
            aln = _OUTDIR + "/2_aligned/palm_extended_{rank}/{taxon}.aln.fasta",
        output:
            trimmed = _OUTDIR + "/3_trimmed/palm_extended_{rank}/{taxon}.trimmed.fasta",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        conda:
            "../envs/trimal.yaml"
        log:
            err = "log/101b_RDRP_phylum_tree/3_trimmed/palm_extended_{rank}_{taxon}.err",
        params:
            resoverlap = config["rdrp_tree"]["trimal_resoverlap"],
            seqoverlap = config["rdrp_tree"]["trimal_seqoverlap"],
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.trimmed}) $(dirname {log.err})
            trimal -in {input.aln} -out {output.trimmed} \
                -resoverlap {params.resoverlap} -seqoverlap {params.seqoverlap} 2> {log.err}
            """


    rule fasttree_ml_palm_extended:
        input:
            trimmed = _OUTDIR + "/3_trimmed/palm_extended_{rank}/{taxon}.trimmed.fasta",
        output:
            nwk = _OUTDIR + "/4_fasttree/palm_extended_{rank}/{taxon}.nwk",
        wildcard_constraints:
            rank = "phylum|class|order|family",
        conda:
            "../envs/fasttree.yaml"
        log:
            log = "log/101b_RDRP_phylum_tree/4_fasttree/palm_extended_{rank}_{taxon}.log",
            err = "log/101b_RDRP_phylum_tree/4_fasttree/palm_extended_{rank}_{taxon}.err",
        params:
            flags = config["rdrp_tree"]["fasttree_flags"],
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        shell:
            """
            mkdir -p $(dirname {output.nwk}) $(dirname {log.log})
            export OMP_NUM_THREADS={threads}
            FastTree {params.flags} {input.trimmed} > {output.nwk} 2> {log.log}
            """
