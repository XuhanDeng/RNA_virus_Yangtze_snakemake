# ESvirtu RdRP identification pipeline
#
# Runs RdRP identification on the ESvirtu database reference sequences that were
# DETECTED in this study's samples (i.e. the non-'acc:' rows of workflow 04's
# merged read_count table — 'acc:' rows are our own assembled contigs and are
# handled by workflow 03).
#
# Step 1: extract the ESvirtu reference rows (first column NOT starting 'acc:')
#         from all_samples.detected_virus.assembly_summary.read_count.tsv.
# Step 2: use those Assembly IDs (read_count col1) to pull the matching metadata
#         rows from virus_pathogen_database.all_metadata.tsv (Assembly = col15).
#         One Assembly may map to multiple Accession rows (multi-segment genomes).
# Step 3: seqkit grep the nucleotide sequences from virus_pathogen_database.fna
#         using the metadata Accession column (col1) as the bait list.
# Step 4+: run the SAME RdRP identification pipeline as workflow 03 on those
#         extracted sequences (RdRpCATCH + ORFfinder + LucaProt -> combine ->
#         motif search -> tier FAAs). Motif profiles are REUSED from workflow 03
#         (result/03_RDRP_identification/4_motif_search/profiles/...), so
#         workflow 03 must have built those profiles first.

configfile: "config/config.yaml"

_ESV_LOG    = "log/04_03_esvirtu_rdrp_identification"
_ESV_OUTDIR = "result/04_03_esvirtu_rdrp_identification"

_MERGE_DIR      = "result/04_modified_esvirtue_ribodetector/es/Merge"
_READCOUNT_TSV  = _MERGE_DIR + "/all_samples.detected_virus.assembly_summary.read_count.tsv"

_ESV_DB_DIR     = "database/esviritu/v3.2.4"
_METADATA_TSV   = _ESV_DB_DIR + "/virus_pathogen_database.all_metadata.tsv"
_DB_FNA         = _ESV_DB_DIR + "/virus_pathogen_database.fna"

_ESV_REF_TSV    = _ESV_OUTDIR + "/1_extract/esvirtu_reference_rows.read_count.tsv"
_ESV_META_TSV   = _ESV_OUTDIR + "/2_metadata/esvirtu_reference.metadata.tsv"
_ESV_FNA        = _ESV_OUTDIR + "/3_sequences/esvirtu_reference.fna"

# ── RdRP identification (mirrors workflow 03) ────────────────────────────────
_STEM           = "esvirtu_reference"
_RC_DIR         = _ESV_OUTDIR + "/4_rdrpcatch"
_ORF_DIR        = _ESV_OUTDIR + "/5_orffinder"
_LP_DIR         = _ESV_OUTDIR + "/6_lucaprot"
_ESV_MOTIF_DIR  = _ESV_OUTDIR + "/7_motif_search"
_ESV_ALL_PROTEINS = _ESV_MOTIF_DIR + "/all_candidates.faa"

# Motif profiles reused from workflow 03 (built there, not rebuilt here).
_MOTIF_OUTDIR   = "result/03_RDRP_identification/4_motif_search"
_MOTIF_SEQ_LIB  = "database/RVMT/RVMT_Zenodo_V4/RVMT_Zenodo_V4/RdRPs/Motifs/Sequence_Library"
_MOTIF_PREFIX   = config["motif_search"]["prefix"]
_MOTIFS         = [1, 2, 3, 4]
_MOTIF_TOOLS    = ["hmmsearch", "psiblast", "mmseqs", "diamond"]

# LucaProt threshold strings (same convention as workflow 03).
_LUCAPROT_THRESHOLD_STR    = "{:.6f}".format(config["lucaprot_rdrp"]["threshold"])
_LUCAPROT_FILTER_THRESHOLD = config["lucaprot_rdrp"].get("filter_threshold", None)
_LUCAPROT_INPUT_STR        = "{:.6f}".format(0.5)
_LUCAPROT_ACTIVE_STR       = (
    "{:.6f}".format(_LUCAPROT_FILTER_THRESHOLD)
    if _LUCAPROT_FILTER_THRESHOLD is not None
    else _LUCAPROT_THRESHOLD_STR
)


rule all:
    input:
        _ESV_REF_TSV,
        _ESV_META_TSV,
        _ESV_FNA,
        # Step 4-5 — RdRpCATCH + LucaProt candidate proteins
        _RC_DIR + "/" + _STEM + "_rdrpcatch_output_annotated.tsv",
        _LP_DIR + "/" + _STEM + "_lucaprot_proteins.faa",
        # Step 6 — combined candidate proteins
        _ESV_ALL_PROTEINS,
        # Step 7 — motif search + tier assignment
        expand(_ESV_MOTIF_DIR + "/search/{tool}/mot.{motif}.tsv",
               tool=_MOTIF_TOOLS, motif=_MOTIFS),
        _ESV_MOTIF_DIR + "/motif_results/tier1a_ABCD.faa",
        _ESV_MOTIF_DIR + "/motif_results/tier1b_CABD_depermuted.faa",
        _ESV_MOTIF_DIR + "/motif_results/tier2_3motif.faa",
        _ESV_MOTIF_DIR + "/motif_results/tier3_no_motif.faa",
        _ESV_MOTIF_DIR + "/motif_results/tier_summary.tsv",


# ── Step 1: extract ESvirtu reference rows ───────────────────────────────────
# Keep only rows whose col4 (Accession) exists in the ESvirtu metadata database
# (col1 of all_metadata.tsv). This is the definitive filter: if col4 matches a
# database accession it is an ESvirtu reference; otherwise it is an assembled contig.

rule extract_esvirtu_reference_rows:
    input:
        tsv      = _READCOUNT_TSV,
        metadata = _METADATA_TSV,
    output:
        tsv = _ESV_REF_TSV,
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        err = _ESV_LOG + "/1_extract/extract_esvirtu_reference_rows.err",
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log.err})
        # Load accessions from metadata col1 into a lookup set.
        # Keep read-count rows (NR>1) where col4 is in that set.
        awk -F'\\t' 'NR==FNR{{acc[$1]=1; next}} FNR>1 && ($4 in acc)' \
            {input.metadata} {input.tsv} > {output.tsv} 2> {log.err}
        echo "ESvirtu reference rows: $(wc -l < {output.tsv})" >> {log.err}
        """


# ── Step 2: pull matching metadata rows (Assembly col1 of step 1 -> col15) ────

rule extract_esvirtu_metadata:
    input:
        ref_tsv  = _ESV_REF_TSV,
        metadata = _METADATA_TSV,
    output:
        tsv = _ESV_META_TSV,
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        err = _ESV_LOG + "/2_metadata/extract_esvirtu_metadata.err",
    shell:
        """
        mkdir -p $(dirname {output.tsv}) $(dirname {log.err})
        # Load Assembly IDs (col1 of step 1) into a set, then keep metadata rows
        # whose Assembly column (col15) is in that set. Header (FNR==1) is dropped.
        awk -F'\\t' 'NR==FNR{{ids[$1]=1; next}} FNR>1 && ($15 in ids)' \
            {input.ref_tsv} {input.metadata} > {output.tsv} 2> {log.err}
        """


# ── Step 3: extract nucleotide sequences by Accession (metadata col1) ─────────

rule extract_esvirtu_sequences:
    input:
        metadata = _ESV_META_TSV,
        fna      = _DB_FNA,
    output:
        fna = _ESV_FNA,
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
        err = _ESV_LOG + "/3_sequences/extract_esvirtu_sequences.err",
    shell:
        """
        mkdir -p $(dirname {output.fna}) $(dirname {log.err})
        # Bait = Accession column (col1) of the matched metadata rows.
        # Prefix all sequence IDs with "esvirtu_" for unambiguous source identification.
        cut -f1 {input.metadata} \
            | seqkit grep --pattern-file /dev/stdin --threads {threads} {input.fna} \
            | seqkit replace --pattern "^" --replacement "esvirtu_" --threads {threads} \
            > {output.fna} 2> {log.err}
        """


# ══════════════════════════════════════════════════════════════════════════════
# Steps 4-7: RdRP identification on the extracted ESvirtu reference sequences
# Mirrors workflow 03 (RdRpCATCH + ORFfinder + LucaProt -> cat -> motif search
# -> tier FAAs). Motif profiles are REUSED from workflow 03's outputs.
# ══════════════════════════════════════════════════════════════════════════════

# ── Step 4: RdRpCATCH on the extracted nucleotide sequences ──────────────────

rule rdrpcatch_esvirtu:
    input:
        _ESV_FNA,
    output:
        tsv        = _RC_DIR + "/" + _STEM + "_rdrpcatch_output_annotated.tsv",
        aa_fasta   = _RC_DIR + "/" + _STEM + "_rdrpcatch_fasta/" + _STEM + "_trimmed_aminoacid_sequences.fasta",
        full_fasta = _RC_DIR + "/" + _STEM + "_rdrpcatch_fasta/" + _STEM + "_full_aminoacid_sequences.fasta",
    log:
        log = _ESV_LOG + "/4_rdrpcatch/rdrpcatch.log",
        err = _ESV_LOG + "/4_rdrpcatch/rdrpcatch.err",
    params:
        output_dir = _RC_DIR,
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
        mkdir -p {params.output_dir} $(dirname {log.log})
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


_ESV_OUTDIR = "result/04_03_esvirtu_rdrp_identification"
_STEM="esvirtu_reference"
_RC_DIR="result/04_03_esvirtu_rdrp_identification" + "/4_rdrpcatch" 
rule extract_rdrpcatch_full_proteins_esvirtu:
    input:
        tsv      = _RC_DIR + "/" + _STEM + "_rdrpcatch_output_annotated.tsv",
        full_faa = _RC_DIR + "/" + _STEM + "_rdrpcatch_fasta/" + _STEM + "_full_aminoacid_sequences.fasta",
    output:
        faa = _RC_DIR + "/" + _STEM + "_rdrpcatch_full_proteins.faa",
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
        err = _ESV_LOG + "/4_rdrpcatch/extract_full.err",
    shell:
        """
        mkdir -p $(dirname {output.faa}) $(dirname {log.err})
        tail -n +2 {input.tsv} \
            | cut -f2 \
            | seqkit grep --pattern-file /dev/stdin --threads {threads} {input.full_faa} \
            > {output.faa} 2> {log.err}
        """


# ── Step 5: ORFfinder + LucaProt ─────────────────────────────────────────────

rule orffinder_esvirtu:
    input:
        _ESV_FNA,
    output:
        aa = _ORF_DIR + "/" + _STEM + "_orfs.faa",
    log:
        log = _ESV_LOG + "/5_orffinder/orffinder.log",
        err = _ESV_LOG + "/5_orffinder/orffinder.err",
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
        out_dir    = _ORF_DIR,
        min_length = config["orffinder"]["min_length"],
        strand     = config["orffinder"]["strand"],
    shell:
        """
        mkdir -p {params.out_dir} $(dirname {log.log})
        ORFfinder -in {input} \
            -ml {params.min_length} \
            -strand {params.strand} \
            -g 1 \
            -s 2 \
            -out {output.aa} \
            -outfmt 0 \
            > {log.log} 2> {log.err}
        """


rule lucaprot_esvirtu:
    input:
        fasta  = _ORF_DIR + "/" + _STEM + "_orfs.faa",
        marker = config["lucaprot"]["marker_db"],
    output:
        csv = _LP_DIR + "/" + _STEM + "_lucaprot_rdrp.csv",
    log:
        log = _ESV_LOG + "/6_lucaprot/lucaprot.log",
        err = _ESV_LOG + "/6_lucaprot/lucaprot.err",
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
        output_dir            = _LP_DIR,
    shell:
        """
        mkdir -p {params.output_dir} $(dirname {log.log})
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

    rule filter_lucaprot_threshold_esvirtu:
        input:
            done = _LP_DIR + "/" + _STEM + "_lucaprot_rdrp.csv",
        output:
            csv = _LP_DIR + "/RdRPs_only_using_threshold" + _LUCAPROT_ACTIVE_STR + ".csv",
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
            rdrp_csv  = _LP_DIR + "/RdRPs_only_using_threshold" + _LUCAPROT_INPUT_STR + ".csv",
            threshold = _LUCAPROT_FILTER_THRESHOLD,
        log:
            out = _ESV_LOG + "/6_lucaprot/filter.log",
            err = _ESV_LOG + "/6_lucaprot/filter.err",
        shell:
            """
            mkdir -p $(dirname {output.csv}) $(dirname {log.out})
            python scripts/03_RDRP_identification/filter_lucaprot_threshold.py \
                --input     {params.rdrp_csv} \
                --output    {output.csv} \
                --threshold {params.threshold} \
                > {log.out} 2> {log.err}
            """


rule extract_lucaprot_proteins_esvirtu:
    input:
        csv = _LP_DIR + "/RdRPs_only_using_threshold" + _LUCAPROT_ACTIVE_STR + ".csv",
    output:
        faa = _LP_DIR + "/" + _STEM + "_lucaprot_proteins.faa",
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
        out = _ESV_LOG + "/6_lucaprot/extract.log",
        err = _ESV_LOG + "/6_lucaprot/extract.err",
    shell:
        """
        mkdir -p $(dirname {output.faa}) $(dirname {log.out})
        python scripts/03_RDRP_identification/extract_lucaprot_proteins_persample.py \
            --csv    {input.csv} \
            --output {output.faa} \
            > {log.out} 2> {log.err}
        """


# ── Step 6: combine RdRpCATCH + LucaProt proteins ────────────────────────────

rule cat_all_candidate_proteins_esvirtu:
    input:
        rc = _RC_DIR + "/" + _STEM + "_rdrpcatch_full_proteins.faa",
        lp = _LP_DIR + "/" + _STEM + "_lucaprot_proteins.faa",
    output:
        faa = _ESV_ALL_PROTEINS,
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    log:
        err = _ESV_LOG + "/7_motif_search/cat_proteins.err",
    shell:
        """
        mkdir -p {_ESV_MOTIF_DIR} $(dirname {log.err})
        cat {input.rc} {input.lp} > {output.faa} 2> {log.err}
        """


# ── Step 7: motif search (reuse workflow 03 profiles) + tier assignment ──────

rule motif_search_hmmsearch_esvirtu:
    input:
        faa = _ESV_ALL_PROTEINS,
        hmm = _MOTIF_OUTDIR + "/profiles/mot.{motif}/HMMfiles/profiles.hmm",
    output:
        tsv = _ESV_MOTIF_DIR + "/search/hmmsearch/mot.{motif}.tsv",
    params:
        raw = _ESV_MOTIF_DIR + "/search/hmmsearch/mot.{motif}_raw.tsv",
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
        log = _ESV_LOG + "/7_motif_search/5_hmmsearch/mot{motif}.log",
        err = _ESV_LOG + "/7_motif_search/5_hmmsearch/mot{motif}.err",
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


rule motif_search_psiblast_esvirtu:
    input:
        faa     = _ESV_ALL_PROTEINS,
        msa_dir = _MOTIF_SEQ_LIB + "/mot.{motif}",
    output:
        tsv = _ESV_MOTIF_DIR + "/search/psiblast/mot.{motif}.tsv",
    params:
        db_dir  = _ESV_MOTIF_DIR + "/search/psiblast/mot.{motif}_blastdb",
        tmp_dir = _ESV_MOTIF_DIR + "/search/psiblast/mot.{motif}_tmp",
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
        log = _ESV_LOG + "/7_motif_search/6_psiblast/mot{motif}.log",
        err = _ESV_LOG + "/7_motif_search/6_psiblast/mot{motif}.err",
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


rule motif_search_mmseqs_esvirtu:
    input:
        faa     = _ESV_ALL_PROTEINS,
        mm_done = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/.build_done",
    output:
        tsv = _ESV_MOTIF_DIR + "/search/mmseqs/mot.{motif}.tsv",
    params:
        query_db           = _ESV_MOTIF_DIR + "/search/mmseqs/mot.{motif}_querydb/DB",
        db_dir             = _MOTIF_OUTDIR + "/profiles/mot.{motif}/MMseqs2_profiles/profiles",
        tmp_dir            = _ESV_MOTIF_DIR + "/search/mmseqs/mot.{motif}_tmp",
        result_dir         = _ESV_MOTIF_DIR + "/search/mmseqs/mot.{motif}_results",
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
        log = _ESV_LOG + "/7_motif_search/7_mmseqs/mot{motif}.log",
        err = _ESV_LOG + "/7_motif_search/7_mmseqs/mot{motif}.err",
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


rule motif_search_diamond_esvirtu:
    input:
        faa     = _ESV_ALL_PROTEINS,
        msa_dir = _MOTIF_SEQ_LIB + "/mot.{motif}",
    output:
        tsv = _ESV_MOTIF_DIR + "/search/diamond/mot.{motif}.tsv",
    params:
        db_dir = _ESV_MOTIF_DIR + "/search/diamond/mot.{motif}_db",
        query  = _ESV_MOTIF_DIR + "/search/diamond/mot.{motif}_query.faa",
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
        log = _ESV_LOG + "/7_motif_search/8_diamond/mot{motif}.log",
        err = _ESV_LOG + "/7_motif_search/8_diamond/mot{motif}.err",
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


rule motif_merge_filter_esvirtu:
    input:
        search_done = expand(_ESV_MOTIF_DIR + "/search/{tool}/mot.{motif}.tsv",
                             tool=_MOTIF_TOOLS, motif=_MOTIFS),
        faa         = _ESV_ALL_PROTEINS,
    output:
        hits       = _ESV_MOTIF_DIR + "/motif_results/motif_hits_best.tsv",
        order      = _ESV_MOTIF_DIR + "/motif_results/motif_order.tsv",
        complete   = _ESV_MOTIF_DIR + "/motif_results/motif_order_complete.tsv",
        depermd    = _ESV_MOTIF_DIR + "/motif_results/sequences_depermuted.faa",
        tier1a     = _ESV_MOTIF_DIR + "/motif_results/tier1a_ABCD.faa",
        tier1b     = _ESV_MOTIF_DIR + "/motif_results/tier1b_CABD_depermuted.faa",
        tier2      = _ESV_MOTIF_DIR + "/motif_results/tier2_3motif.faa",
        tier3      = _ESV_MOTIF_DIR + "/motif_results/tier3_no_motif.faa",
        suspicious = _ESV_MOTIF_DIR + "/motif_results/motif_order_suspicious.tsv",
        summary    = _ESV_MOTIF_DIR + "/motif_results/tier_summary.tsv",
    params:
        search_dir = _ESV_MOTIF_DIR + "/search",
        outdir     = _ESV_MOTIF_DIR + "/motif_results",
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
        log = _ESV_LOG + "/7_motif_search/9_merge/merge.log",
        err = _ESV_LOG + "/7_motif_search/9_merge/merge.err",
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
