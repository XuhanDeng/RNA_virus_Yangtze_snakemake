# RdRP identification pipeline
# Step 1: RdRpCATCH — HMM-based RdRP candidate detection (high sensitivity, uses all databases)
# Step 2: ORFfinder + LucaProt — deep learning RdRp detection on all ORFs
# Step 3: Concatenate RdRpCATCH + LucaProt proteins into one FAA
# Step 4: RdRp motif A-D search (RVMT Sequence_Library) — hmmsearch, psiblast, diamond

configfile: "config/config.yaml"

MIN_LEN  = config["seqkit"]["min_length"]
ICTV_STEM = "Riboviria_sequences"

RDRP_LOG_DIR  = "log/03_RDRP_identification"
RDRP_DIR      = "result/03_RDRP_identification/1_RdRpCATCH"
ORF_DIR       = "result/03_RDRP_identification/2_orffinder"
LP_DIR        = "result/03_RDRP_identification/3_lucaprot"

# ── ICTV path constants ────────────────────────────────────────────────────────
_ICTV_RC_DIR      = RDRP_DIR + "/ICTV"
_ICTV_LP_DIR      = LP_DIR + "/ICTV"
_ICTV_ORF_DIR     = ORF_DIR + "/ICTV"
_ICTV_MOTIF_DIR   = "result/03_RDRP_identification/4_motif_search/ICTV"
_ICTV_ALL_PROTEINS = _ICTV_MOTIF_DIR + "/all_candidates.faa"


def ictv_targets():
    if not config["ICTV"]["use"]:
        return []
    return [
        _ICTV_RC_DIR + "/" + ICTV_STEM + "_rdrpcatch_output_annotated.tsv",
        _ICTV_LP_DIR + "/ICTV_lucaprot_rdrp.csv",
        _ICTV_LP_DIR + "/ICTV_lucaprot_proteins.faa",
        expand(_ICTV_MOTIF_DIR + "/search/{tool}/mot.{motif}.tsv",
               tool=_MOTIF_TOOLS, motif=_MOTIFS),
        _ICTV_MOTIF_DIR + "/motif_results/motif_hits_best.tsv",
        _ICTV_MOTIF_DIR + "/motif_results/motif_order.tsv",
        _ICTV_MOTIF_DIR + "/motif_results/motif_order_complete.tsv",
        _ICTV_MOTIF_DIR + "/motif_results/sequences_depermuted.faa",
        _ICTV_MOTIF_DIR + "/motif_results/tier1a_ABCD.faa",
        _ICTV_MOTIF_DIR + "/motif_results/tier1b_CABD_depermuted.faa",
        _ICTV_MOTIF_DIR + "/motif_results/tier2_3motif.faa",
        _ICTV_MOTIF_DIR + "/motif_results/tier3_no_motif.faa",
        _ICTV_MOTIF_DIR + "/motif_results/motif_order_suspicious.tsv",
        _ICTV_MOTIF_DIR + "/motif_results/tier_summary.tsv",
    ]


_LUCAPROT_THRESHOLD_STR     = "{:.6f}".format(config["lucaprot_rdrp"]["threshold"])
_LUCAPROT_FILTER_THRESHOLD  = config["lucaprot_rdrp"].get("filter_threshold", None)
_LUCAPROT_INPUT_STR         = "{:.6f}".format(0.5)
_LUCAPROT_ACTIVE_STR        = (
    "{:.6f}".format(_LUCAPROT_FILTER_THRESHOLD)
    if _LUCAPROT_FILTER_THRESHOLD is not None
    else _LUCAPROT_THRESHOLD_STR
)

_MOTIF_OUTDIR  = "result/03_RDRP_identification/4_motif_search"
_MOTIF_SEQ_LIB = "database/RVMT/RVMT_Zenodo_V4/RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library"
_MOTIF_PREFIX  = config["motif_search"]["prefix"]
_MOTIFS        = [1, 2, 3, 4]
_MOTIF_TOOLS   = ["hmmsearch", "psiblast", "mmseqs", "diamond"]
_ALL_PROTEINS  = _MOTIF_OUTDIR + "/all_candidates.faa"


rule all:
    input:
        expand(
            RDRP_DIR + "/{sample}/{sample}_scaffolds_rename_" + str(MIN_LEN) + "_rdrpcatch_output_annotated.tsv",
            sample=config["rna_samples"]
        ),
        ictv_targets(),
        expand(
            LP_DIR + "/{sample}/{sample}_lucaprot_proteins.faa",
            sample=config["rna_samples"]
        ),
        # Step 3 — combined protein FAA
        _ALL_PROTEINS,
        # Step 4 — motif search results
        expand(_MOTIF_OUTDIR + "/profiles/mot.{motif}/HMMfiles/profiles.hmm",   motif=_MOTIFS),
        expand(_MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/.build_done", motif=_MOTIFS),
        expand(_MOTIF_OUTDIR + "/search/{tool}/mot.{motif}.tsv",
               tool=_MOTIF_TOOLS, motif=_MOTIFS),
        # Step 5 — motif merge + filter
        _MOTIF_OUTDIR + "/motif_results/motif_hits_best.tsv",
        _MOTIF_OUTDIR + "/motif_results/motif_order.tsv",
        _MOTIF_OUTDIR + "/motif_results/motif_order_complete.tsv",
        _MOTIF_OUTDIR + "/motif_results/sequences_depermuted.faa",
        _MOTIF_OUTDIR + "/motif_results/tier1a_ABCD.faa",
        _MOTIF_OUTDIR + "/motif_results/tier1b_CABD_depermuted.faa",
        _MOTIF_OUTDIR + "/motif_results/tier2_3motif.faa",
        _MOTIF_OUTDIR + "/motif_results/tier3_no_motif.faa",
        _MOTIF_OUTDIR + "/motif_results/motif_order_suspicious.tsv",
        _MOTIF_OUTDIR + "/motif_results/tier_summary.tsv",


# ── Step 1: RdRpCATCH (samples) ───────────────────────────────────────────────

rule rdrpcatch:#ok
    input:
        config["rna_reformated_scaffolds_dir"] + "/rename_" + str(MIN_LEN) + "/{sample}_scaffolds_rename_" + str(MIN_LEN) + ".fasta"
    output:
        tsv       = RDRP_DIR + "/{sample}/{sample}_scaffolds_rename_" + str(MIN_LEN) + "_rdrpcatch_output_annotated.tsv",
        aa_fasta  = RDRP_DIR + "/{sample}/{sample}_scaffolds_rename_" + str(MIN_LEN) + "_rdrpcatch_fasta/{sample}_scaffolds_rename_" + str(MIN_LEN) + "_trimmed_aminoacid_sequences.fasta",
        full_fasta = RDRP_DIR + "/{sample}/{sample}_scaffolds_rename_" + str(MIN_LEN) + "_rdrpcatch_fasta/{sample}_scaffolds_rename_" + str(MIN_LEN) + "_full_aminoacid_sequences.fasta",
    log:
        log = RDRP_LOG_DIR + "/rdrpcatch/{sample}.log",
        err = RDRP_LOG_DIR + "/rdrpcatch/{sample}.err"
    params:
        output_dir = RDRP_DIR + "/{sample}",
        seq_type   = config["rdrp_catch"]["seq_type"],
        db_dir     = config["rdrp_catch"]["db_dir"],
        db_options = config["rdrp_catch"]["db_options"]
    threads: config["rdrp_catch"]["threads"]
    resources:
        mem_mb_per_cpu  = config["regular_memory"],
        runtime         = config["rdrp_catch"]["runtime"],
        cpus_per_task   = config["rdrp_catch"]["threads"],
        slurm_partition = config["regular_partition"],
        slurm_account   = config["account"]
    conda:
        "../envs/rdrp_catch.yaml"
    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p $(dirname {log.log})
        rdrpcatch scan --input {input} \
                  --output {params.output_dir} \
                  --cpus {threads} \
                  --seq-type {params.seq_type} \
                  --db-dir {params.db_dir} \
                  --db-options {params.db_options} \
                  --extended-output \
                  --overwrite \
                  > {log.log} 2> {log.err}
        """


# ── Step 2: ORFfinder + LucaProt (samples) ───────────────────────────────────

rule orffinder:#ok
    input:
        config["rna_reformated_scaffolds_dir"] + "/rename_" + str(MIN_LEN) + "/{sample}_scaffolds_rename_" + str(MIN_LEN) + ".fasta"
    output:
        aa = ORF_DIR + "/{sample}/{sample}_orfs.faa"
    log:
        log = RDRP_LOG_DIR + "/orffinder/{sample}.log",
        err = RDRP_LOG_DIR + "/orffinder/{sample}.err"
    conda:
        "../envs/orffinder.yaml"
    threads: config["orffinder"]["threads"]
    resources:
        mem_mb_per_cpu  = config["orffinder"]["memory"],
        runtime         = config["orffinder"]["runtime"],
        cpus_per_task   = config["orffinder"]["threads"],
        slurm_partition = config["orffinder"]["partition"],
        slurm_account   = config["orffinder"]["account"]
    params:
        out_dir    = ORF_DIR + "/{sample}",
        min_length = config["orffinder"]["min_length"],
        strand     = config["orffinder"]["strand"]
    shell:
        """
        mkdir -p {params.out_dir}
        mkdir -p $(dirname {log.log})
        ORFfinder -in {input} \
            -ml {params.min_length} \
            -strand {params.strand} \
            -g 1 \
            -s 2 \
            -out {output.aa} \
            -outfmt 0 \
            > {log.log} 2> {log.err}
        """

rule lucaprot_rdrp:#ok
    input:
        fasta  = ORF_DIR + "/{sample}/{sample}_orfs.faa",
        marker = config["lucaprot"]["marker_db"]
    output:
        csv = LP_DIR + "/{sample}/{sample}_lucaprot_rdrp.csv"
    log:
        log = RDRP_LOG_DIR + "/lucaprot/{sample}.log",
        err = RDRP_LOG_DIR + "/lucaprot/{sample}.err"
    conda:
        "../envs/lucaprot.yaml"
    threads: config["lucaprot_rdrp"]["threads"]
    resources:
        mem_mb_per_cpu  = config["lucaprot_rdrp"]["memory"],
        runtime         = config["lucaprot_rdrp"]["runtime"],
        cpus_per_task   = config["lucaprot_rdrp"]["threads"],
        slurm_partition = config["lucaprot_rdrp"]["partition"],
        slurm_account   = config["lucaprot_rdrp"]["account"],
        slurm_extra     = "'--gpus-per-task={}'".format(config["lucaprot_rdrp"]["gpus"]) if config["lucaprot_rdrp"]["gpu_id"] >= 0 else ""
    params:
        install_dir           = config["lucaprot"]["install_dir"],
        db_dir                = config["lucaprot"]["db_dir"],
        gpu_id                = config["lucaprot_rdrp"]["gpu_id"],
        gpu_devices           = config["lucaprot_rdrp"]["gpu_devices"],
        truncation_seq_length = config["lucaprot_rdrp"]["truncation_seq_length"],
        dataset_name          = "rdrp_40_extend",
        dataset_type          = "protein",
        task_type             = "binary_class",
        model_type            = "sefn",
        time_str              = config["lucaprot_rdrp"]["time_str"],
        step                  = config["lucaprot_rdrp"]["step"],
        threshold             = config["lucaprot_rdrp"]["threshold"],
        print_per_number      = config["lucaprot_rdrp"]["print_per_number"],
        output_dir            = LP_DIR + "/{sample}"
    shell:
        """
        mkdir -p {params.output_dir}
        mkdir -p $(dirname {log.log})
        fasta_abs=$(realpath {input.fasta})
        out_abs=$(realpath {output.csv})
        db_abs=$(realpath {params.db_dir})
        log_abs=$(realpath {log.log})
        err_abs=$(realpath {log.err})

        if [ "{params.gpu_id}" -ge 0 ]; then
            export CUDA_VISIBLE_DEVICES="{params.gpu_devices}"
        fi

        cd {params.install_dir}/src
        python predict_many_samples.py \
            --fasta_file $fasta_abs \
            --save_file $out_abs \
            --truncation_seq_length {params.truncation_seq_length} \
            --dataset_name {params.dataset_name} \
            --dataset_type {params.dataset_type} \
            --task_type {params.task_type} \
            --model_type {params.model_type} \
            --time_str {params.time_str} \
            --step {params.step} \
            --threshold {params.threshold} \
            --print_per_number {params.print_per_number} \
            --gpu_id {params.gpu_id} \
            --torch_hub_dir $db_abs \
            > $log_abs 2> $log_abs
        """

# Step 2b: filter lucaprot CSV to higher threshold (optional)
if _LUCAPROT_FILTER_THRESHOLD is not None:

    rule filter_lucaprot_threshold:
        wildcard_constraints:
            sample = r"[^/]+"
        input:
            csv = LP_DIR + "/{sample}/RdRPs_only_using_threshold" + _LUCAPROT_INPUT_STR + ".csv",
        output:
            csv = LP_DIR + "/{sample}/RdRPs_only_using_threshold" + _LUCAPROT_ACTIVE_STR + ".csv",
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
            threshold = _LUCAPROT_FILTER_THRESHOLD,
        log:
            out = "log/03_RDRP_identification/lucaprot_filter/{sample}.log",
            err = "log/03_RDRP_identification/lucaprot_filter/{sample}.err",
        shell:
            """
            mkdir -p $(dirname {output.csv}) log/03_RDRP_identification/lucaprot_filter
            python scripts/03_RDRP_identification/filter_lucaprot_threshold.py \
                --input     {input.csv} \
                --output    {output.csv} \
                --threshold {params.threshold} \
                > {log.out} 2> {log.err}
            """

# Step 2c: extract LucaProt proteins per sample
rule extract_lucaprot_proteins_persample:#ok
    input:
        csv = LP_DIR + "/{sample}/RdRPs_only_using_threshold" + _LUCAPROT_ACTIVE_STR + ".csv",
    output:
        faa = LP_DIR + "/{sample}/{sample}_lucaprot_proteins.faa",
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = 1,
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        out = "log/03_RDRP_identification/lucaprot_extract/{sample}.log",
        err = "log/03_RDRP_identification/lucaprot_extract/{sample}.err",
    shell:
        """
        mkdir -p $(dirname {output.faa}) log/03_RDRP_identification/lucaprot_extract
        python scripts/03_RDRP_identification/extract_lucaprot_proteins_persample.py \
            --csv    {input.csv} \
            --output {output.faa} \
            > {log.out} 2> {log.err}
        """


# ── Step 3: Extract full-length RdRpCATCH proteins + concatenate with LucaProt ─

rule extract_rdrpcatch_full_proteins:#ok
    input:
        tsv      = RDRP_DIR + "/{sample}/{sample}_scaffolds_rename_" + str(MIN_LEN) + "_rdrpcatch_output_annotated.tsv",
        full_faa = RDRP_DIR + "/{sample}/{sample}_scaffolds_rename_" + str(MIN_LEN) + "_rdrpcatch_fasta/{sample}_scaffolds_rename_" + str(MIN_LEN) + "_full_aminoacid_sequences.fasta",
    output:
        faa = RDRP_DIR + "/{sample}/{sample}_rdrpcatch_full_proteins.faa",
    conda:
        "../envs/seqkit-spade.yaml"
    threads: config["seqkit"]["threads"]
    resources:
        mem_mb_per_cpu  = config["seqkit"]["memory"],
        runtime         = config["seqkit"]["runtime"],
        cpus_per_task   = config["seqkit"]["threads"],
        slurm_partition = config["seqkit"]["partition"],
        slurm_account   = config["seqkit"]["account"],
    log:
        err = "log/03_RDRP_identification/extract_rc_full/{sample}.err",
    shell:
        """
        mkdir -p $(dirname {output.faa}) log/03_RDRP_identification/extract_rc_full
        tail -n +2 {input.tsv} \
            | cut -f2 \
            | seqkit grep --pattern-file /dev/stdin --threads {threads} {input.full_faa} \
            > {output.faa} 2> {log.err}
        """


rule cat_all_candidate_proteins:#ok
    input:
        rc = expand(
            RDRP_DIR + "/{sample}/{sample}_rdrpcatch_full_proteins.faa",
            sample=config["rna_samples"]
        ),
        lp = expand(
            LP_DIR + "/{sample}/{sample}_lucaprot_proteins.faa",
            sample=config["rna_samples"]
        ),
    output:
        faa = _ALL_PROTEINS,
    conda:
        "../envs/python.yaml"
    threads: 1
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = 1,
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        err = "log/03_RDRP_identification/motif_search/cat_proteins.err",
    shell:
        """
        mkdir -p {_MOTIF_OUTDIR} log/03_RDRP_identification/motif_search
        cat {input.rc} {input.lp} > {output.faa} 2> {log.err}
        """


# ── Step 4: RdRp motif A-D search (RVMT Sequence_Library) ────────────────────

# ── HMM (hmmer) ───────────────────────────────────────────────────────────────

rule motif_add_consensus:#ok
    input:
        msa_dir = _MOTIF_SEQ_LIB + "/mot.{motif}",
    output:
        done = _MOTIF_OUTDIR + "/profiles/mot.{motif}/msaFiles/.consensus_done",
    params:
        out_msa = _MOTIF_OUTDIR + "/profiles/mot.{motif}/msaFiles",
    conda:
        "../envs/hhsuite.yaml"
    threads: config["motif_search"]["threads"]
    resources:
        mem_mb_per_cpu  = config["motif_search"]["memory"],
        runtime         = config["motif_search"]["runtime"],
        cpus_per_task   = config["motif_search"]["threads"],
        slurm_partition = config["motif_search"]["partition"],
        slurm_account   = config["motif_search"]["account"],
    log:
        log = RDRP_LOG_DIR + "/4_motif_search/1_consensus/mot{motif}.log",
        err = RDRP_LOG_DIR + "/4_motif_search/1_consensus/mot{motif}.err",
    shell:
        """
        module load {config[motif_search][parallel_module]}
        mkdir -p {params.out_msa}
        cp {input.msa_dir}/*.afa {params.out_msa}/
        parallel -j {threads} \
            hhconsensus -M 50 -cov 50 \
                -i {{}} \
                -oa2m {params.out_msa}/{{/.}}.Cons.msa.afa \
            ::: {params.out_msa}/*_woCon.afa \
        >> {log.log} 2>> {log.err}
        parallel -j {threads} sed -i '1d' {{}} \
            ::: {params.out_msa}/*.Cons.msa.afa \
        >> {log.log} 2>> {log.err}
        touch {output.done}
        """


rule motif_build_hmm:#ok,check
    input:
        msa_dir = _MOTIF_SEQ_LIB + "/mot.{motif}",
    output:
        hmm_db = _MOTIF_OUTDIR + "/profiles/mot.{motif}/HMMfiles/profiles.hmm",
    params:
        hmm_dir         = _MOTIF_OUTDIR + "/profiles/mot.{motif}/HMMfiles",
        log_dir         = _MOTIF_OUTDIR + "/profiles/mot.{motif}/logs",
        parallel_module = config["motif_search"]["parallel_module"],
    conda:
        "../envs/hmmer.yaml"
    threads: config["motif_search"]["threads"]
    resources:
        mem_mb_per_cpu  = config["motif_search"]["memory"],
        runtime         = config["motif_search"]["runtime"],
        cpus_per_task   = config["motif_search"]["threads"],
        slurm_partition = config["motif_search"]["partition"],
        slurm_account   = config["motif_search"]["account"],
    log:
        log = RDRP_LOG_DIR + "/4_motif_search/2_hmmbuild/mot{motif}.log",
        err = RDRP_LOG_DIR + "/4_motif_search/2_hmmbuild/mot{motif}.err",
    shell:
        """
        mkdir -p {params.hmm_dir} {params.log_dir}
        if module load {params.parallel_module} 2>/dev/null && command -v parallel &>/dev/null; then
            echo "Using GNU parallel ({threads} jobs)" >> {log.log}
            parallel -j {threads} \
                hmmbuild --informat afa \
                    -n {{/.}} \
                    -o {params.log_dir}/{{/.}}.hmm.log \
                    {params.hmm_dir}/{{/.}}.hmm \
                    {{}} \
                ::: {input.msa_dir}/*.afa \
            >> {log.log} 2>> {log.err}
        else
            echo "parallel not available, running sequentially" >> {log.log}
            for f in {input.msa_dir}/*.afa; do
                name=$(basename "$f" .afa)
                hmmbuild --informat afa \
                    -n "$name" \
                    -o {params.log_dir}/"$name".hmm.log \
                    {params.hmm_dir}/"$name".hmm \
                    "$f" \
                >> {log.log} 2>> {log.err}
            done
        fi
        cat {params.hmm_dir}/*.hmm > {output.hmm_db}
        """


rule motif_search_hmmsearch:#ok,check
    input:
        faa = _ALL_PROTEINS,
        hmm = _MOTIF_OUTDIR + "/profiles/mot.{motif}/HMMfiles/profiles.hmm",
    output:
        tsv = _MOTIF_OUTDIR + "/search/hmmsearch/mot.{motif}.tsv",
    params:
        raw = _MOTIF_OUTDIR + "/search/hmmsearch/mot.{motif}_raw.tsv",
    conda:
        "../envs/hmmer.yaml"
    threads: config["motif_search"]["threads"]
    resources:
        mem_mb_per_cpu  = config["motif_search"]["memory"],
        runtime         = config["motif_search"]["runtime"],
        cpus_per_task   = config["motif_search"]["threads"],
        slurm_partition = config["motif_search"]["partition"],
        slurm_account   = config["motif_search"]["account"],
    log:
        log = RDRP_LOG_DIR + "/4_motif_search/5_hmmsearch/mot{motif}.log",
        err = RDRP_LOG_DIR + "/4_motif_search/5_hmmsearch/mot{motif}.err",
    shell:
        """
        mkdir -p {_MOTIF_OUTDIR}/search/hmmsearch
        hmmsearch \
            --noali \
            --cpu {threads} \
            -E 0.5 \
            --incE 0.5 \
            --domtblout {params.raw} \
            {input.hmm} \
            {input.faa} \
        >> {log.log} 2>> {log.err}
        grep -v "^#" {params.raw} > {output.tsv} 2>> {log.err} || true
        """


# ── MMseqs2 ───────────────────────────────────────────────────────────────────

rule motif_msa_to_stockholm:#ok
    input:
        msa_dir = _MOTIF_SEQ_LIB + "/mot.{motif}",
    output:
        sto_done = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/sto/.sto_done",
    params:
        sto_dir = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/sto",
    conda:
        "../envs/mafft_hmmer_seqkit.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        log = RDRP_LOG_DIR + "/4_motif_search/3a_sto/mot{motif}.log",
        err = RDRP_LOG_DIR + "/4_motif_search/3a_sto/mot{motif}.err",
    shell:
        """
        for f in {input.msa_dir}/*.afa; do
            name=$(basename "$f" .afa)
            esl-reformat stockholm "$f" \
                > {params.sto_dir}/"$name".sto \
                2>> {log.err}
        done
        touch {output.sto_done}
        """


rule motif_build_mmseqs_profile:#ok
    input:
        sto_done = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/sto/.sto_done",
    output:
        mm_done = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/.build_done",
    params:
        sto_dir = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/sto",
        db_dir  = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/profiles",
        parallel_module = config["motif_search"]["parallel_module"],
    conda:
        "../envs/mmseqs2.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        log = RDRP_LOG_DIR + "/4_motif_search/3b_mmseqs_profile/mot{motif}.log",
        err = RDRP_LOG_DIR + "/4_motif_search/3b_mmseqs_profile/mot{motif}.err",
    shell:
        """
        mkdir -p {params.db_dir}
        build_one() {{
            f="$1"
            db_dir="$2"
            log="$3"
            err="$4"
            name=$(basename "$f" .sto)
            mmseqs convertmsa "$f" \
                "$db_dir"/"$name"_msa_db \
            >> "$log" 2>> "$err"
            mmseqs msa2profile \
                "$db_dir"/"$name"_msa_db \
                "$db_dir"/"$name"_profile \
                --match-mode 1 \
            >> "$log" 2>> "$err"
        }}
        export -f build_one
        if module load {params.parallel_module} 2>/dev/null && command -v parallel &>/dev/null; then
            echo "Using GNU parallel ({threads} jobs)" >> {log.log}
            parallel -j {threads} build_one {{}} {params.db_dir} {log.log} {log.err} \
                ::: {params.sto_dir}/*.sto
        else
            echo "parallel not available, running sequentially" >> {log.log}
            for f in {params.sto_dir}/*.sto; do
                build_one "$f" {params.db_dir} {log.log} {log.err}
            done
        fi
        touch {output.mm_done}
        """


rule motif_search_mmseqs:#ok
    input:
        faa     = _ALL_PROTEINS,
        mm_done = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/.build_done",
    output:
        tsv = _MOTIF_OUTDIR + "/search/mmseqs/mot.{motif}.tsv",
    params:
        query_db        = _MOTIF_OUTDIR + "/search/mmseqs/mot.{motif}_querydb/DB",
        db_dir          = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/profiles",
        tmp_dir         = _MOTIF_OUTDIR + "/search/mmseqs/mot.{motif}_tmp",
        result_dir      = _MOTIF_OUTDIR + "/search/mmseqs/mot.{motif}_results",
        parallel_module    = config["motif_search"]["parallel_module"],
        split_memory_limit = config["motif_search"]["split_memory_limit"],
        par_jobs           = 4,
    conda:
        "../envs/mmseqs2.yaml"
    threads: config["motif_search"]["threads"]
    resources:
        mem_mb_per_cpu  = config["motif_search"]["memory"],
        runtime         = config["motif_search"]["runtime"],
        cpus_per_task   = config["motif_search"]["threads"],
        slurm_partition = config["motif_search"]["partition"],
        slurm_account   = config["motif_search"]["account"],
    log:
        log = RDRP_LOG_DIR + "/4_motif_search/7_mmseqs/mot{motif}.log",
        err = RDRP_LOG_DIR + "/4_motif_search/7_mmseqs/mot{motif}.err",
    shell:
        """
        mkdir -p $(dirname {params.query_db}) {params.tmp_dir} {params.result_dir}
        mmseqs createdb {input.faa} {params.query_db} >> {log.log} 2>> {log.err}
        search_one() {{
            profile="$1"
            query_db="$2"
            result_dir="$3"
            tmp_dir="$4"
            log="$5"
            err="$6"
            threads_per_job="$7"
            name=$(basename "$profile" _profile)
            result="$result_dir"/"$name"_result/DB
            mkdir -p $(dirname "$result")
            mmseqs search \
                "$profile" \
                "$query_db" \
                "$result" \
                "$tmp_dir"/"$name"_tmp \
                -s 7.5 -e 0.5 \
                --threads "$threads_per_job" \
                --split-memory-limit {params.split_memory_limit} \
                --max-seqs 6000 \
            >> "$log" 2>> "$err"
            mmseqs convertalis \
                "$profile" \
                "$query_db" \
                "$result" \
                "$result_dir"/"$name".tsv \
                --format-output \
                    "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,qframe,mismatch,qcov,tcov" \
            >> "$log" 2>> "$err"
        }}
        export -f search_one
        THREADS_PER_JOB=$(( {threads} / {params.par_jobs} ))
        [ "$THREADS_PER_JOB" -lt 1 ] && THREADS_PER_JOB=1
        if module load {params.parallel_module} 2>/dev/null && command -v parallel &>/dev/null; then
            echo "Using GNU parallel ({params.par_jobs} jobs, $THREADS_PER_JOB threads each)" >> {log.log}
            parallel -j {params.par_jobs} search_one {{}} \
                {params.query_db} {params.result_dir} {params.tmp_dir} \
                {log.log} {log.err} "$THREADS_PER_JOB" \
                ::: {params.db_dir}/*_profile
        else
            echo "parallel not available, running sequentially" >> {log.log}
            for profile in {params.db_dir}/*_profile; do
                search_one "$profile" \
                    {params.query_db} {params.result_dir} {params.tmp_dir} \
                    {log.log} {log.err} {threads}
            done
        fi
        cat {params.result_dir}/*.tsv > {output.tsv} 2>> {log.err}
        rm -rf {params.result_dir}/*_result
        """


# ── PSI-BLAST + Diamond ───────────────────────────────────────────────────────

rule motif_search_psiblast:
    input:
        faa     = _ALL_PROTEINS,
        msa_dir = _MOTIF_SEQ_LIB + "/mot.{motif}",
    output:
        tsv = _MOTIF_OUTDIR + "/search/psiblast/mot.{motif}.tsv",
    params:
        db_dir  = _MOTIF_OUTDIR + "/search/psiblast/mot.{motif}_blastdb",
        tmp_dir = _MOTIF_OUTDIR + "/search/psiblast/mot.{motif}_tmp",
        stem    = _MOTIF_PREFIX,
    conda:
        "../envs/blast.yaml"
    threads: config["motif_search"]["threads"]
    resources:
        mem_mb_per_cpu  = config["motif_search"]["memory"],
        runtime         = config["motif_search"]["runtime"],
        cpus_per_task   = config["motif_search"]["threads"],
        slurm_partition = config["motif_search"]["partition"],
        slurm_account   = config["motif_search"]["account"],
    log:
        log = RDRP_LOG_DIR + "/4_motif_search/6_psiblast/mot{motif}.log",
        err = RDRP_LOG_DIR + "/4_motif_search/6_psiblast/mot{motif}.err",
    shell:
        """
        mkdir -p {params.db_dir} {params.tmp_dir}
        makeblastdb \
            -in {input.faa} \
            -dbtype prot \
            -title {params.stem} \
            -out {params.db_dir}/{params.stem} \
        >> {log.log} 2>> {log.err}
        for cluster in {input.msa_dir}/*_woCon.afa; do
            cluster_name=$(basename "$cluster" "_woCon.afa")
            n=$(grep -c "^>" "$cluster" 2>/dev/null || echo 0)
            if [ "$n" -lt 2 ]; then
                echo "Skipping $cluster_name: only $n sequence(s)" >> {log.log}
                continue
            fi
            psiblast \
                -word_size 2 \
                -evalue 0.5 \
                -max_target_seqs 100000 \
                -threshold 9 \
                -dbsize 20000000 \
                -in_msa "$cluster" \
                -qcov_hsp_perc 1 \
                -num_threads {threads} \
                -db {params.db_dir}/{params.stem} \
                -outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore" \
                -out {params.tmp_dir}/"${{cluster_name}}"_out.tsv \
            >> {log.log} 2>> {log.err}
            awk -v cls="$cluster_name" \
                'BEGIN{{OFS="\t"}} {{print $0, cls}}' \
                {params.tmp_dir}/"${{cluster_name}}"_out.tsv \
            >> {output.tsv}
        done
        """


rule motif_search_diamond:
    input:
        faa     = _ALL_PROTEINS,
        msa_dir = _MOTIF_SEQ_LIB + "/mot.{motif}",
    output:
        tsv = _MOTIF_OUTDIR + "/search/diamond/mot.{motif}.tsv",
    params:
        db_dir   = _MOTIF_OUTDIR + "/search/diamond/mot.{motif}_db",
        query    = _MOTIF_OUTDIR + "/search/diamond/mot.{motif}_query.faa",
        stem     = _MOTIF_PREFIX,
    conda:
        "../envs/diamond.yaml"
    threads: config["motif_search"]["threads"]
    resources:
        mem_mb_per_cpu  = config["motif_search"]["memory"],
        runtime         = config["motif_search"]["runtime"],
        cpus_per_task   = config["motif_search"]["threads"],
        slurm_partition = config["motif_search"]["partition"],
        slurm_account   = config["motif_search"]["account"],
    log:
        log = RDRP_LOG_DIR + "/4_motif_search/8_diamond/mot{motif}.log",
        err = RDRP_LOG_DIR + "/4_motif_search/8_diamond/mot{motif}.err",
    shell:
        """
        mkdir -p {params.db_dir} $(dirname {params.query})
        cat {input.msa_dir}/*_woCon.afa \
            | grep -v "^-" \
            | awk '/^>/{{header=$0}} !/^>/{{gsub(/-/,""); print header; print}}' \
            > {params.query} 2>> {log.err}
        diamond makedb \
            --in {input.faa} \
            -d {params.db_dir}/{params.stem} \
            -p {threads} \
        >> {log.log} 2>> {log.err}
        diamond blastp \
            -q {params.query} \
            --db {params.db_dir}/{params.stem} \
            -p {threads} \
            -b 3.0 \
            -k 0 \
            --query-cover 50 \
            --id 60 \
            --more-sensitive \
            --evalue 0.05 \
            --outfmt 6 qseqid sseqid evalue gapopen pident length \
                        qstart qend sstart send qlen slen mismatch bitscore \
            -o {output.tsv} \
        >> {log.log} 2>> {log.err}
        """


# ── Step 5: Merge + filter motif hits ────────────────────────────────────────

rule motif_merge_filter:
    input:
        search_done = expand(_MOTIF_OUTDIR + "/search/{tool}/mot.{motif}.tsv",
                             tool=_MOTIF_TOOLS, motif=_MOTIFS),
        faa         = _ALL_PROTEINS,
    output:
        hits       = _MOTIF_OUTDIR + "/motif_results/motif_hits_best.tsv",
        order      = _MOTIF_OUTDIR + "/motif_results/motif_order.tsv",
        complete   = _MOTIF_OUTDIR + "/motif_results/motif_order_complete.tsv",
        depermd    = _MOTIF_OUTDIR + "/motif_results/sequences_depermuted.faa",
        tier1a     = _MOTIF_OUTDIR + "/motif_results/tier1a_ABCD.faa",
        tier1b     = _MOTIF_OUTDIR + "/motif_results/tier1b_CABD_depermuted.faa",
        tier2      = _MOTIF_OUTDIR + "/motif_results/tier2_3motif.faa",
        tier3      = _MOTIF_OUTDIR + "/motif_results/tier3_no_motif.faa",
        suspicious = _MOTIF_OUTDIR + "/motif_results/motif_order_suspicious.tsv",
        summary    = _MOTIF_OUTDIR + "/motif_results/tier_summary.tsv",
    params:
        search_dir = _MOTIF_OUTDIR + "/search",
        outdir     = _MOTIF_OUTDIR + "/motif_results",
        tools      = " ".join(_MOTIF_TOOLS),
        motifs     = " ".join(str(m) for m in _MOTIFS),
    conda:
        "../envs/pandas_biopython.yaml"
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        log = RDRP_LOG_DIR + "/4_motif_search/9_merge/merge.log",
        err = RDRP_LOG_DIR + "/4_motif_search/9_merge/merge.err",
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log.log})
        python scripts/03_RDRP_identification/motif_merge_filter.py \
            --search_dir {params.search_dir} \
            --input_faa  {input.faa} \
            --outdir     {params.outdir} \
            --tools      {params.tools} \
            --motifs     {params.motifs} \
        > {log.log} 2> {log.err}
        """


# ── ICTV: Full RdRP pipeline for ICTV Riboviria reference sequences ───────────

if config["ICTV"]["use"]:

    rule rdrpcatch_ICTV:
        input:
            config["ICTV"]["fasta"]
        output:
            tsv       = _ICTV_RC_DIR + "/" + ICTV_STEM + "_rdrpcatch_output_annotated.tsv",
            aa_fasta  = _ICTV_RC_DIR + "/" + ICTV_STEM + "_rdrpcatch_fasta/" + ICTV_STEM + "_trimmed_aminoacid_sequences.fasta",
            full_fasta = _ICTV_RC_DIR + "/" + ICTV_STEM + "_rdrpcatch_fasta/" + ICTV_STEM + "_full_aminoacid_sequences.fasta",
        log:
            log = "log/03_RDRP_identification/rdrpcatch/ICTV.log",
            err = "log/03_RDRP_identification/rdrpcatch/ICTV.err",
        params:
            output_dir = _ICTV_RC_DIR,
            seq_type   = config["rdrp_catch"]["seq_type"],
            db_dir     = config["rdrp_catch"]["db_dir"],
            db_options = config["rdrp_catch"]["db_options"],
        threads: config["rdrp_catch"]["threads"]
        resources:
            mem_mb_per_cpu  = config["regular_memory"],
            runtime         = config["rdrp_catch"]["runtime"],
            cpus_per_task   = config["rdrp_catch"]["threads"],
            slurm_partition = config["regular_partition"],
            slurm_account   = config["account"],
        conda:
            "../envs/rdrp_catch.yaml"
        shell:
            """
            mkdir -p {params.output_dir}
            mkdir -p $(dirname {log.log})
            rdrpcatch scan --input {input} \
                      --output {params.output_dir} \
                      --cpus {threads} \
                      --seq-type {params.seq_type} \
                      --db-dir {params.db_dir} \
                      --db-options {params.db_options} \
                      --extended-output \
                      --overwrite \
                      > {log.log} 2> {log.err}
            """

    rule orffinder_ICTV:
        input:
            config["ICTV"]["fasta"]
        output:
            aa = _ICTV_ORF_DIR + "/ICTV_orfs.faa",
        log:
            log = "log/03_RDRP_identification/orffinder/ICTV.log",
            err = "log/03_RDRP_identification/orffinder/ICTV.err",
        conda:
            "../envs/orffinder.yaml"
        threads: config["orffinder"]["threads"]
        resources:
            mem_mb_per_cpu  = config["orffinder"]["memory"],
            runtime         = config["orffinder"]["runtime"],
            cpus_per_task   = config["orffinder"]["threads"],
            slurm_partition = config["orffinder"]["partition"],
            slurm_account   = config["orffinder"]["account"],
        params:
            out_dir    = _ICTV_ORF_DIR,
            min_length = config["orffinder"]["min_length"],
            strand     = config["orffinder"]["strand"],
        shell:
            """
            mkdir -p {params.out_dir}
            mkdir -p $(dirname {log.log})
            ORFfinder -in {input} \
                -ml {params.min_length} \
                -strand {params.strand} \
                -g 1 \
                -s 2 \
                -out {output.aa} \
                -outfmt 0 \
                > {log.log} 2> {log.err}
            """

    rule lucaprot_ICTV:
        input:
            fasta  = _ICTV_ORF_DIR + "/ICTV_orfs.faa",
            marker = config["lucaprot"]["marker_db"],
        output:
            csv = _ICTV_LP_DIR + "/ICTV_lucaprot_rdrp.csv",
        log:
            log = "log/03_RDRP_identification/lucaprot/ICTV.log",
            err = "log/03_RDRP_identification/lucaprot/ICTV.err",
        conda:
            "../envs/lucaprot.yaml"
        threads: config["lucaprot_rdrp"]["threads"]
        resources:
            mem_mb_per_cpu  = config["lucaprot_rdrp"]["memory"],
            runtime         = config["lucaprot_rdrp"]["runtime"],
            cpus_per_task   = config["lucaprot_rdrp"]["threads"],
            slurm_partition = config["lucaprot_rdrp"]["partition"],
            slurm_account   = config["lucaprot_rdrp"]["account"],
            slurm_extra     = "'--gpus-per-task={}'".format(config["lucaprot_rdrp"]["gpus"]) if config["lucaprot_rdrp"]["gpu_id"] >= 0 else "",
        params:
            install_dir           = config["lucaprot"]["install_dir"],
            db_dir                = config["lucaprot"]["db_dir"],
            gpu_id                = config["lucaprot_rdrp"]["gpu_id"],
            gpu_devices           = config["lucaprot_rdrp"]["gpu_devices"],
            truncation_seq_length = config["lucaprot_rdrp"]["truncation_seq_length"],
            dataset_name          = "rdrp_40_extend",
            dataset_type          = "protein",
            task_type             = "binary_class",
            model_type            = "sefn",
            time_str              = config["lucaprot_rdrp"]["time_str"],
            step                  = config["lucaprot_rdrp"]["step"],
            threshold             = config["lucaprot_rdrp"]["threshold"],
            print_per_number      = config["lucaprot_rdrp"]["print_per_number"],
            output_dir            = _ICTV_LP_DIR,
        shell:
            """
            mkdir -p {params.output_dir}
            mkdir -p $(dirname {log.log})
            fasta_abs=$(realpath {input.fasta})
            out_abs=$(realpath {output.csv})
            db_abs=$(realpath {params.db_dir})
            log_abs=$(realpath {log.log})
            err_abs=$(realpath {log.err})

            if [ "{params.gpu_id}" -ge 0 ]; then
                export CUDA_VISIBLE_DEVICES="{params.gpu_devices}"
            fi

            cd {params.install_dir}/src
            python predict_many_samples.py \
                --fasta_file $fasta_abs \
                --save_file $out_abs \
                --truncation_seq_length {params.truncation_seq_length} \
                --dataset_name {params.dataset_name} \
                --dataset_type {params.dataset_type} \
                --task_type {params.task_type} \
                --model_type {params.model_type} \
                --time_str {params.time_str} \
                --step {params.step} \
                --threshold {params.threshold} \
                --print_per_number {params.print_per_number} \
                --gpu_id {params.gpu_id} \
                --torch_hub_dir $db_abs \
                > $log_abs 2> $err_abs
            """

    if _LUCAPROT_FILTER_THRESHOLD is not None:

        rule filter_lucaprot_threshold_ICTV:
            input:
                csv = _ICTV_LP_DIR + "/RdRPs_only_using_threshold" + _LUCAPROT_INPUT_STR + ".csv",
            output:
                csv = _ICTV_LP_DIR + "/RdRPs_only_using_threshold" + _LUCAPROT_ACTIVE_STR + ".csv",
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
                threshold = _LUCAPROT_FILTER_THRESHOLD,
            log:
                out = "log/03_RDRP_identification/lucaprot_filter/ICTV.log",
                err = "log/03_RDRP_identification/lucaprot_filter/ICTV.err",
            shell:
                """
                mkdir -p $(dirname {output.csv}) log/03_RDRP_identification/lucaprot_filter
                python scripts/03_RDRP_identification/filter_lucaprot_threshold.py \
                    --input     {input.csv} \
                    --output    {output.csv} \
                    --threshold {params.threshold} \
                    > {log.out} 2> {log.err}
                """

    rule extract_lucaprot_proteins_ICTV:
        input:
            csv = _ICTV_LP_DIR + "/RdRPs_only_using_threshold" + _LUCAPROT_ACTIVE_STR + ".csv",
        output:
            faa = _ICTV_LP_DIR + "/ICTV_lucaprot_proteins.faa",
        conda:
            "../envs/python.yaml"
        threads: 1
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = 1,
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        log:
            out = "log/03_RDRP_identification/lucaprot_extract/ICTV.log",
            err = "log/03_RDRP_identification/lucaprot_extract/ICTV.err",
        shell:
            """
            mkdir -p $(dirname {output.faa}) log/03_RDRP_identification/lucaprot_extract
            python scripts/03_RDRP_identification/extract_lucaprot_proteins_persample.py \
                --csv    {input.csv} \
                --output {output.faa} \
                > {log.out} 2> {log.err}
            """

    rule extract_rdrpcatch_full_proteins_ICTV:
        input:
            tsv      = _ICTV_RC_DIR + "/" + ICTV_STEM + "_rdrpcatch_output_annotated.tsv",
            full_faa = _ICTV_RC_DIR + "/" + ICTV_STEM + "_rdrpcatch_fasta/" + ICTV_STEM + "_full_aminoacid_sequences.fasta",
        output:
            faa = _ICTV_RC_DIR + "/ICTV_rdrpcatch_full_proteins.faa",
        conda:
            "../envs/seqkit-spade.yaml"
        threads: config["seqkit"]["threads"]
        resources:
            mem_mb_per_cpu  = config["seqkit"]["memory"],
            runtime         = config["seqkit"]["runtime"],
            cpus_per_task   = config["seqkit"]["threads"],
            slurm_partition = config["seqkit"]["partition"],
            slurm_account   = config["seqkit"]["account"],
        log:
            err = "log/03_RDRP_identification/extract_rc_full/ICTV.err",
        shell:
            """
            mkdir -p $(dirname {output.faa}) log/03_RDRP_identification/extract_rc_full
            tail -n +2 {input.tsv} \
                | cut -f2 \
                | seqkit grep --pattern-file /dev/stdin --threads {threads} {input.full_faa} \
                > {output.faa} 2> {log.err}
            """

    rule cat_all_candidate_proteins_ICTV:
        input:
            rc = _ICTV_RC_DIR + "/ICTV_rdrpcatch_full_proteins.faa",
            lp = _ICTV_LP_DIR + "/ICTV_lucaprot_proteins.faa",
        output:
            faa = _ICTV_ALL_PROTEINS,
        threads: 1
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = 1,
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        log:
            err = "log/03_RDRP_identification/motif_search/ICTV_cat_proteins.err",
        shell:
            """
            mkdir -p {_ICTV_MOTIF_DIR} log/03_RDRP_identification/motif_search
            cat {input.rc} {input.lp} > {output.faa} 2> {log.err}
            """

    # ── ICTV motif search — reuse profiles built from _MOTIF_SEQ_LIB ─────────

    rule motif_search_hmmsearch_ICTV:
        input:
            faa = _ICTV_ALL_PROTEINS,
            hmm = _MOTIF_OUTDIR + "/profiles/mot.{motif}/HMMfiles/profiles.hmm",
        output:
            tsv = _ICTV_MOTIF_DIR + "/search/hmmsearch/mot.{motif}.tsv",
        params:
            raw = _ICTV_MOTIF_DIR + "/search/hmmsearch/mot.{motif}_raw.tsv",
        conda:
            "../envs/hmmer.yaml"
        threads: config["motif_search"]["threads"]
        resources:
            mem_mb_per_cpu  = config["motif_search"]["memory"],
            runtime         = config["motif_search"]["runtime"],
            cpus_per_task   = config["motif_search"]["threads"],
            slurm_partition = config["motif_search"]["partition"],
            slurm_account   = config["motif_search"]["account"],
        log:
            log = RDRP_LOG_DIR + "/4_motif_search/5_hmmsearch/ICTV_mot{motif}.log",
            err = RDRP_LOG_DIR + "/4_motif_search/5_hmmsearch/ICTV_mot{motif}.err",
        shell:
            """
            mkdir -p $(dirname {output.tsv})
            hmmsearch \
                --tblout {params.raw} \
                --cpu {threads} \
                -E 0.5 \
                {input.hmm} {input.faa} \
            >> {log.log} 2>> {log.err}
            grep -v "^#" {params.raw} \
                | awk '{{print $1, $3, $5, $6}}' OFS='\t' \
                > {output.tsv}
            """

    rule motif_search_psiblast_ICTV:
        input:
            faa     = _ICTV_ALL_PROTEINS,
            msa_dir = _MOTIF_SEQ_LIB + "/mot.{motif}",
        output:
            tsv = _ICTV_MOTIF_DIR + "/search/psiblast/mot.{motif}.tsv",
        params:
            db_dir  = _ICTV_MOTIF_DIR + "/search/psiblast/mot.{motif}_blastdb",
            tmp_dir = _ICTV_MOTIF_DIR + "/search/psiblast/mot.{motif}_tmp",
            stem    = _MOTIF_PREFIX,
        conda:
            "../envs/blast.yaml"
        threads: config["motif_search"]["threads"]
        resources:
            mem_mb_per_cpu  = config["motif_search"]["memory"],
            runtime         = config["motif_search"]["runtime"],
            cpus_per_task   = config["motif_search"]["threads"],
            slurm_partition = config["motif_search"]["partition"],
            slurm_account   = config["motif_search"]["account"],
        log:
            log = RDRP_LOG_DIR + "/4_motif_search/6_psiblast/ICTV_mot{motif}.log",
            err = RDRP_LOG_DIR + "/4_motif_search/6_psiblast/ICTV_mot{motif}.err",
        shell:
            """
            mkdir -p {params.db_dir} {params.tmp_dir}
            makeblastdb \
                -in {input.faa} \
                -dbtype prot \
                -title {params.stem} \
                -out {params.db_dir}/{params.stem} \
            >> {log.log} 2>> {log.err}
            for cluster in {input.msa_dir}/*_woCon.afa; do
                cluster_name=$(basename "$cluster" "_woCon.afa")
                n=$(grep -c "^>" "$cluster" 2>/dev/null || echo 0)
                if [ "$n" -lt 2 ]; then
                    echo "Skipping $cluster_name: only $n sequence(s)" >> {log.log}
                    continue
                fi
                psiblast \
                    -word_size 2 \
                    -evalue 0.5 \
                    -max_target_seqs 100000 \
                    -threshold 9 \
                    -dbsize 20000000 \
                    -in_msa "$cluster" \
                    -qcov_hsp_perc 1 \
                    -num_threads {threads} \
                    -db {params.db_dir}/{params.stem} \
                    -outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore" \
                    -out {params.tmp_dir}/"${{cluster_name}}"_out.tsv \
                >> {log.log} 2>> {log.err}
                awk -v cls="$cluster_name" \
                    'BEGIN{{OFS="\t"}} {{print $0, cls}}' \
                    {params.tmp_dir}/"${{cluster_name}}"_out.tsv \
                >> {output.tsv}
            done
            """

    rule motif_search_mmseqs_ICTV:
        input:
            faa     = _ICTV_ALL_PROTEINS,
            mm_done = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/.build_done",
        output:
            tsv = _ICTV_MOTIF_DIR + "/search/mmseqs/mot.{motif}.tsv",
        params:
            query_db           = _ICTV_MOTIF_DIR + "/search/mmseqs/mot.{motif}_querydb/DB",
            db_dir             = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/profiles",
            tmp_dir            = _ICTV_MOTIF_DIR + "/search/mmseqs/mot.{motif}_tmp",
            result_dir         = _ICTV_MOTIF_DIR + "/search/mmseqs/mot.{motif}_results",
            parallel_module    = config["motif_search"]["parallel_module"],
            split_memory_limit = config["motif_search"]["split_memory_limit"],
            par_jobs           = 4,
        conda:
            "../envs/mmseqs2.yaml"
        threads: config["motif_search"]["threads"]
        resources:
            mem_mb_per_cpu  = config["motif_search"]["memory"],
            runtime         = config["motif_search"]["runtime"],
            cpus_per_task   = config["motif_search"]["threads"],
            slurm_partition = config["motif_search"]["partition"],
            slurm_account   = config["motif_search"]["account"],
        log:
            log = RDRP_LOG_DIR + "/4_motif_search/7_mmseqs/ICTV_mot{motif}.log",
            err = RDRP_LOG_DIR + "/4_motif_search/7_mmseqs/ICTV_mot{motif}.err",
        shell:
            """
            mkdir -p $(dirname {params.query_db}) {params.tmp_dir} {params.result_dir}
            mmseqs createdb {input.faa} {params.query_db} >> {log.log} 2>> {log.err}
            search_one() {{
                profile="$1"
                query_db="$2"
                result_dir="$3"
                tmp_dir="$4"
                log="$5"
                err="$6"
                threads_per_job="$7"
                name=$(basename "$profile" _profile)
                result="$result_dir"/"$name"_result/DB
                mkdir -p $(dirname "$result")
                mmseqs search \
                    "$profile" \
                    "$query_db" \
                    "$result" \
                    "$tmp_dir"/"$name"_tmp \
                    -s 7.5 -e 0.5 \
                    --threads "$threads_per_job" \
                    --split-memory-limit {params.split_memory_limit} \
                    --max-seqs 6000 \
                >> "$log" 2>> "$err"
                mmseqs convertalis \
                    "$profile" \
                    "$query_db" \
                    "$result" \
                    "$result_dir"/"$name".tsv \
                    --format-output \
                        "query,target,evalue,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen" \
                >> "$log" 2>> "$err"
            }}
            export -f search_one
            THREADS_PER_JOB=$(( {threads} / {params.par_jobs} ))
            [ "$THREADS_PER_JOB" -lt 1 ] && THREADS_PER_JOB=1
            if module load {params.parallel_module} 2>/dev/null && command -v parallel &>/dev/null; then
                echo "Using GNU parallel ({params.par_jobs} jobs, $THREADS_PER_JOB threads each)" >> {log.log}
                parallel -j {params.par_jobs} search_one {{}} \
                    {params.query_db} {params.result_dir} {params.tmp_dir} \
                    {log.log} {log.err} "$THREADS_PER_JOB" \
                    ::: {params.db_dir}/*_profile
            else
                echo "parallel not available, running sequentially" >> {log.log}
                for profile in {params.db_dir}/*_profile; do
                    search_one "$profile" \
                        {params.query_db} {params.result_dir} {params.tmp_dir} \
                        {log.log} {log.err} {threads}
                done
            fi
            cat {params.result_dir}/*.tsv > {output.tsv} 2>> {log.err}
            rm -rf {params.result_dir}/*_result
            """

    rule motif_search_diamond_ICTV:
        input:
            faa     = _ICTV_ALL_PROTEINS,
            msa_dir = _MOTIF_SEQ_LIB + "/mot.{motif}",
        output:
            tsv = _ICTV_MOTIF_DIR + "/search/diamond/mot.{motif}.tsv",
        params:
            db_dir = _ICTV_MOTIF_DIR + "/search/diamond/mot.{motif}_db",
            query  = _ICTV_MOTIF_DIR + "/search/diamond/mot.{motif}_query.faa",
            stem   = _MOTIF_PREFIX,
        conda:
            "../envs/diamond.yaml"
        threads: config["motif_search"]["threads"]
        resources:
            mem_mb_per_cpu  = config["motif_search"]["memory"],
            runtime         = config["motif_search"]["runtime"],
            cpus_per_task   = config["motif_search"]["threads"],
            slurm_partition = config["motif_search"]["partition"],
            slurm_account   = config["motif_search"]["account"],
        log:
            log = RDRP_LOG_DIR + "/4_motif_search/8_diamond/ICTV_mot{motif}.log",
            err = RDRP_LOG_DIR + "/4_motif_search/8_diamond/ICTV_mot{motif}.err",
        shell:
            """
            mkdir -p {params.db_dir} $(dirname {params.query})
            cat {input.msa_dir}/*_woCon.afa \
                | awk '/^>/{{header=$0}} !/^>/{{gsub(/-/,""); print header; print}}' \
                > {params.query} 2>> {log.err}
            diamond makedb \
                --in {input.faa} \
                -d {params.db_dir}/{params.stem} \
                -p {threads} \
            >> {log.log} 2>> {log.err}
            diamond blastp \
                -q {params.query} \
                --db {params.db_dir}/{params.stem} \
                -p {threads} \
                -b 3.0 \
                -k 0 \
                --query-cover 50 \
                --id 60 \
                --more-sensitive \
                --evalue 0.05 \
                --outfmt 6 qseqid sseqid evalue gapopen pident length \
                            qstart qend sstart send qlen slen mismatch bitscore \
                -o {output.tsv} \
            >> {log.log} 2>> {log.err}
            """

    rule motif_merge_filter_ICTV:
        input:
            search_done = expand(_ICTV_MOTIF_DIR + "/search/{tool}/mot.{motif}.tsv",
                                 tool=_MOTIF_TOOLS, motif=_MOTIFS),
            faa         = _ICTV_ALL_PROTEINS,
        output:
            hits       = _ICTV_MOTIF_DIR + "/motif_results/motif_hits_best.tsv",
            order      = _ICTV_MOTIF_DIR + "/motif_results/motif_order.tsv",
            complete   = _ICTV_MOTIF_DIR + "/motif_results/motif_order_complete.tsv",
            depermd    = _ICTV_MOTIF_DIR + "/motif_results/sequences_depermuted.faa",
            tier1a     = _ICTV_MOTIF_DIR + "/motif_results/tier1a_ABCD.faa",
            tier1b     = _ICTV_MOTIF_DIR + "/motif_results/tier1b_CABD_depermuted.faa",
            tier2      = _ICTV_MOTIF_DIR + "/motif_results/tier2_3motif.faa",
            tier3      = _ICTV_MOTIF_DIR + "/motif_results/tier3_no_motif.faa",
            suspicious = _ICTV_MOTIF_DIR + "/motif_results/motif_order_suspicious.tsv",
            summary    = _ICTV_MOTIF_DIR + "/motif_results/tier_summary.tsv",
        params:
            search_dir = _ICTV_MOTIF_DIR + "/search",
            outdir     = _ICTV_MOTIF_DIR + "/motif_results",
            tools      = " ".join(_MOTIF_TOOLS),
            motifs     = " ".join(str(m) for m in _MOTIFS),
        conda:
            "../envs/pandas_biopython.yaml"
        threads: config["small_job"]["threads"]
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = config["small_job"]["threads"],
            slurm_partition = config["small_job"]["partition"],
            slurm_account   = config["small_job"]["account"],
        log:
            log = RDRP_LOG_DIR + "/4_motif_search/9_merge/ICTV_merge.log",
            err = RDRP_LOG_DIR + "/4_motif_search/9_merge/ICTV_merge.err",
        shell:
            """
            mkdir -p {params.outdir} $(dirname {log.log})
            python scripts/03_RDRP_identification/motif_merge_filter.py \
                --search_dir {params.search_dir} \
                --input_faa  {input.faa} \
                --outdir     {params.outdir} \
                --tools      {params.tools} \
                --motifs     {params.motifs} \
            > {log.log} 2> {log.err}
            """

