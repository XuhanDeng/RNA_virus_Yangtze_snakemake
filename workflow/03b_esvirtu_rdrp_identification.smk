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
#          extracted sequences:
#          RdRpCATCH -> palm_annot -> ORFfinder -> LucaProt -> palm_annot ->
#          merge_rdrp_results -> extract proteins per source (RdRpCATCH / LucaProt).
#          Mirrors the simplified no-tier design from workflow 03 exactly.

configfile: "config/config.yaml"

_ESV_LOG    = "log/03b_esvirtu_rdrp_identification"
_ESV_OUTDIR = "result/03b_esvirtu_rdrp_identification"

_MERGE_DIR     = "result/04_modified_esvirtue_ribodetector/es/Merge"
_READCOUNT_TSV = _MERGE_DIR + "/all_samples.detected_virus.assembly_summary.read_count.tsv"

_ESV_DB_DIR  = "database/esviritu/v3.2.4"
_METADATA_TSV = _ESV_DB_DIR + "/virus_pathogen_database.all_metadata.tsv"
_DB_FNA       = _ESV_DB_DIR + "/virus_pathogen_database.fna"

_ESV_REF_TSV = _ESV_OUTDIR + "/1_extract/esvirtu_reference_rows.read_count.tsv"
_ESV_META_TSV = _ESV_OUTDIR + "/2_metadata/esvirtu_reference.metadata.tsv"
_ESV_FNA      = _ESV_OUTDIR + "/3_sequences/esvirtu_reference.fna"

# ── RdRP identification (mirrors workflow 03) ────────────────────────────────
_STEM        = "esvirtu_reference"
_RC_DIR      = _ESV_OUTDIR + "/4_rdrpcatch"
_ORF_DIR     = _ESV_OUTDIR + "/5_orffinder"
_LP_DIR      = _ESV_OUTDIR + "/6_lucaprot"
_PS_DIR      = _ESV_OUTDIR + "/7_palmscan"
_LPS_DIR     = _ESV_OUTDIR + "/8_lucaprot_palmscan"
_MERGE_OUT   = _ESV_OUTDIR + "/9_merged"
_PROTEIN_DIR = _ESV_OUTDIR + "/10_proteins"
_PALM_DIR    = _ESV_OUTDIR + "/10b_palm_extracted"
_FINAL_DIR   = _ESV_OUTDIR + "/11_final"
_TAXDIR      = _ESV_OUTDIR + "/12_taxonomy"
_PHYLUM_DIR  = _ESV_OUTDIR + "/13_phylum_cluster"

# LucaProt threshold strings (same convention as workflow 03).
_LUCAPROT_THRESHOLD_STR    = "{:.6f}".format(config["lucaprot_rdrp"]["threshold"])
_LUCAPROT_FILTER_THRESHOLD = config["lucaprot_rdrp"].get("filter_threshold", None)
_LUCAPROT_INPUT_STR        = "{:.6f}".format(0.5)
_LUCAPROT_ACTIVE_STR       = (
    "{:.6f}".format(_LUCAPROT_FILTER_THRESHOLD)
    if _LUCAPROT_FILTER_THRESHOLD is not None
    else _LUCAPROT_THRESHOLD_STR
)

# Protein source key (same as workflow 03)
_PROTEIN_SOURCE_KEY = {"RdRpCATCH": "ps", "LucaProt": "lpt"}
_PROTEIN_SOURCE_FASTA = {
    "ps":  _PS_DIR  + "/" + _STEM + "_rdrp_trimmed.faa",
    "lpt": _LPS_DIR + "/" + _STEM + "_lucaprot_rdrp_trimmed.faa",
}


rule all:
    input:
        _ESV_REF_TSV,
        _ESV_META_TSV,
        _ESV_FNA,
        _RC_DIR  + "/" + _STEM + "_rdrpcatch_output_annotated.tsv",
        _PS_DIR  + "/" + _STEM + "_palmscan_hits.tsv",
        _LP_DIR  + "/" + _STEM + "_lucaprot_rdrp.csv",
        _LPS_DIR + "/" + _STEM + "_lucaprot_palmscan_hits.tsv",
        _MERGE_OUT + "/" + _STEM + "_rdrp_merged.tsv",
        _PROTEIN_DIR + "/RdRpCATCH.faa",
        _PROTEIN_DIR + "/LucaProt.faa",
        _PROTEIN_DIR + "/final_proteins_summary.tsv",
        _PALM_DIR + "/RdRpCATCH/palm_regions.tsv",
        _PALM_DIR + "/LucaProt/palm_regions.tsv",
        _FINAL_DIR + "/combined_full.faa",
        _FINAL_DIR + "/combined_palm_core.faa",
        _FINAL_DIR + "/combined_palm_extended.faa",
        _FINAL_DIR + "/combined_palm_regions.tsv",
        expand(_TAXDIR + "/{query}_taxonomy.tsv", query=["full_length", "palm_core", "palm_extended"]),
        expand(_PHYLUM_DIR + "/{query}_phylum", query=["full_length", "palm_core", "palm_extended"]),
        expand(_PHYLUM_DIR + "/full_length_{rank}", rank=["class", "order", "family"]),
        expand(_PHYLUM_DIR + "/{palmset}_{rank}",
               palmset=["palm_core", "palm_extended"], rank=["class", "order", "family"]),


# ── Step 1: extract ESvirtu reference rows ───────────────────────────────────

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
        awk -F'\\t' 'NR==FNR{{acc[$1]=1; next}} FNR>1 && ($4 in acc)' \
            {input.metadata} {input.tsv} > {output.tsv} 2> {log.err}
        echo "ESvirtu reference rows: $(wc -l < {output.tsv})" >> {log.err}
        """


# ── Step 2: pull matching metadata rows ──────────────────────────────────────

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
        awk -F'\\t' 'NR==FNR{{ids[$1]=1; next}} FNR>1 && ($15 in ids)' \
            {input.ref_tsv} {input.metadata} > {output.tsv} 2> {log.err}
        """


# ── Step 3: extract nucleotide sequences ─────────────────────────────────────

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
        cut -f1 {input.metadata} \
            | seqkit grep --pattern-file /dev/stdin --threads {threads} {input.fna} \
            | seqkit replace --pattern "^" --replacement "esvirtu_" --threads {threads} \
            > {output.fna} 2> {log.err}
        """


# ── Step 4: RdRpCATCH ────────────────────────────────────────────────────────

rule rdrpcatch_esvirtu:
    input:
        _ESV_FNA,
    output:
        tsv      = _RC_DIR + "/" + _STEM + "_rdrpcatch_output_annotated.tsv",
        aa_fasta = _RC_DIR + "/" + _STEM + "_rdrpcatch_fasta/" + _STEM + "_trimmed_aminoacid_sequences.fasta",
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


# ── Step 4b: palm_annot on RdRpCATCH trimmed AA ──────────────────────────────

rule palm_annot_esvirtu:
    input:
        _RC_DIR + "/" + _STEM + "_rdrpcatch_fasta/" + _STEM + "_trimmed_aminoacid_sequences.fasta",
    output:
        tsv     = _PS_DIR + "/" + _STEM + "_palmscan_hits.tsv",
        fev     = _PS_DIR + "/" + _STEM + "_palm_annot.fev",
        rdrp_aa = _PS_DIR + "/" + _STEM + "_rdrp_trimmed.faa",
    log:
        log = _ESV_LOG + "/7_palmscan/palm_annot.log",
        err = _ESV_LOG + "/7_palmscan/palm_annot.err",
    params:
        output_dir    = _PS_DIR,
        palm_annot_py = config["palm_annot"]["install_dir"] + "/py/palm_annot.py",
        fev2tsv_py    = config["palm_annot"]["install_dir"] + "/py/fev2tsv.py",
        seqtype       = config["palm_annot"]["seqtype"],
        minscore      = config["palm_annot"]["minscore"],
        minpssmscore  = config["palm_annot"]["minpssmscore"],
        tmpdir        = config["palm_annot"]["tmpdir"] + "/esvirtu",
    threads: config["palm_annot"]["threads"]
    resources:
        mem_mb_per_cpu  = config["palm_annot"]["memory"],
        runtime         = config["palm_annot"]["runtime"],
        cpus_per_task   = config["palm_annot"]["threads"],
        slurm_partition = config["palm_annot"]["partition"],
        slurm_account   = config["palm_annot"]["account"],
    conda:
        "../envs/palm_annot.yaml"
    shell:
        """
        mkdir -p {params.output_dir} {params.tmpdir} $(dirname {log.log})
        python {params.palm_annot_py} \
            --input        {input} \
            --seqtype      {params.seqtype} \
            --fev          {output.fev} \
            --rdrp         {output.rdrp_aa} \
            --minscore     {params.minscore} \
            --minpssmscore {params.minpssmscore} \
            --threads      {threads} \
            --tmpdir       {params.tmpdir} \
            > {log.log} 2> {log.err}
        python {params.fev2tsv_py} \
            --input  {output.fev} \
            --output {output.tsv} \
            --header yes \
            >> {log.log} 2>> {log.err}
        """


# ── Step 5: ORFfinder ────────────────────────────────────────────────────────

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


# ── Step 6: LucaProt ─────────────────────────────────────────────────────────

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
        threads: 1
        resources:
            mem_mb_per_cpu  = config["small_job"]["memory"],
            runtime         = config["small_job"]["runtime"],
            cpus_per_task   = 1,
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
    threads: 1
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = 1,
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


# ── Step 6b: palm_annot on LucaProt proteins ─────────────────────────────────

rule palm_annot_lucaprot_esvirtu:
    input:
        aa = _LP_DIR + "/" + _STEM + "_lucaprot_proteins.faa",
    output:
        tsv     = _LPS_DIR + "/" + _STEM + "_lucaprot_palmscan_hits.tsv",
        fev     = _LPS_DIR + "/" + _STEM + "_lucaprot_palm_annot.fev",
        rdrp_aa = _LPS_DIR + "/" + _STEM + "_lucaprot_rdrp_trimmed.faa",
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
        seqtype       = config["palm_annot"]["seqtype"],
        minscore      = config["palm_annot"]["minscore"],
        minpssmscore  = config["palm_annot"]["minpssmscore"],
        tmpdir        = config["palm_annot"]["tmpdir"] + "/esvirtu_lucaprot",
    log:
        out = _ESV_LOG + "/8_lucaprot_palmscan/palm_annot.log",
        err = _ESV_LOG + "/8_lucaprot_palmscan/palm_annot.err",
    shell:
        """
        mkdir -p {params.tmpdir} $(dirname {log.out})
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


# ── Step 7: merge RdRpCATCH + LucaProt results ───────────────────────────────

rule merge_rdrp_results_esvirtu:
    input:
        rdrpcatch         = _RC_DIR  + "/" + _STEM + "_rdrpcatch_output_annotated.tsv",
        palmscan          = _PS_DIR  + "/" + _STEM + "_palmscan_hits.tsv",
        lucaprot          = _LP_DIR  + "/RdRPs_only_using_threshold" + _LUCAPROT_ACTIVE_STR + ".csv",
        lucaprot_palmscan = _LPS_DIR + "/" + _STEM + "_lucaprot_palmscan_hits.tsv",
    output:
        merged  = _MERGE_OUT + "/" + _STEM + "_rdrp_merged.tsv",
        regions = _MERGE_OUT + "/" + _STEM + "_rdrp_regions.tsv",
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
        script = config["scripts"]["merge_rdrp_results"],
    log:
        out = _ESV_LOG + "/9_merged/merge.log",
        err = _ESV_LOG + "/9_merged/merge.err",
    shell:
        """
        mkdir -p {_MERGE_OUT} $(dirname {log.out})
        python {params.script} \
            --rdrpcatch          {input.rdrpcatch} \
            --palmscan           {input.palmscan} \
            --lucaprot           {input.lucaprot} \
            --lucaprot-palmscan  {input.lucaprot_palmscan} \
            --output             {output.merged} \
            --output-regions     {output.regions} \
            > {log.out} 2> {log.err}
        """


# ── Step 8: write per-source protein ID lists ─────────────────────────────────

rule write_protein_id_lists_esvirtu:
    input:
        merged = _MERGE_OUT + "/" + _STEM + "_rdrp_merged.tsv",
    output:
        expand(_PROTEIN_DIR + "/{cat}_ids.txt", cat=["RdRpCATCH", "LucaProt"]),
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
        script = config["scripts"]["write_protein_id_lists"],
        outdir = _PROTEIN_DIR,
    log:
        out = _ESV_LOG + "/10_proteins/write_ids.log",
        err = _ESV_LOG + "/10_proteins/write_ids.err",
    shell:
        """
        mkdir -p {params.outdir} $(dirname {log.out})
        python {params.script} \
            --merged {input.merged} \
            --outdir {params.outdir} \
            > {log.out} 2> {log.err}
        """


# ── Step 8b: human-readable index of the final protein FASTAs ────────────────

rule summarize_final_proteins_esvirtu:
    input:
        merged = _MERGE_OUT + "/" + _STEM + "_rdrp_merged.tsv",
    output:
        _PROTEIN_DIR + "/final_proteins_summary.tsv",
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
        out = _ESV_LOG + "/10_proteins/final_proteins_summary.log",
        err = _ESV_LOG + "/10_proteins/final_proteins_summary.err",
    shell:
        """
        mkdir -p {_PROTEIN_DIR} $(dirname {log.out})
        python scripts/03_RDRP_identification/summarize_final_proteins.py \
            --merged {input.merged} \
            --output {output} \
            > {log.out} 2> {log.err}
        """


# ── Step 9: extract per-source protein FASTAs ────────────────────────────────

rule extract_proteins_esvirtu:
    input:
        fasta   = lambda wc: _PROTEIN_SOURCE_FASTA[_PROTEIN_SOURCE_KEY[wc.pcat]],
        id_list = _PROTEIN_DIR + "/{pcat}_ids.txt",
    output:
        fasta = _PROTEIN_DIR + "/{pcat}.faa",
    wildcard_constraints:
        pcat = "RdRpCATCH|LucaProt",
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
        err = _ESV_LOG + "/10_proteins/{pcat}.err",
    shell:
        """
        mkdir -p {_PROTEIN_DIR} $(dirname {log.err})
        seqkit grep \
            --pattern-file {input.id_list} \
            --threads {threads} \
            {input.fasta} \
            > {output.fasta} \
            2> {log.err}
        """


rule extract_palm_regions_esvirtu:
    input:
        rdrpcatch_motif = _PROTEIN_DIR + "/RdRpCATCH.faa",
        lucaprot_motif  = _PROTEIN_DIR + "/LucaProt.faa",
        rdrpcatch_seqs  = _RC_DIR + "/" + _STEM + "_rdrpcatch_fasta/" + _STEM + "_trimmed_aminoacid_sequences.fasta",
        lucaprot_seqs   = _ORF_DIR + "/" + _STEM + "_orfs.faa",
        script          = "scripts/03_RDRP_identification/extract_palm_regions.py",
    output:
        rc_tsv      = _PALM_DIR + "/RdRpCATCH/palm_regions.tsv",
        lp_tsv      = _PALM_DIR + "/LucaProt/palm_regions.tsv",
        rc_full     = _PALM_DIR + "/RdRpCATCH/rdrp_full.faa",
        lp_full     = _PALM_DIR + "/LucaProt/rdrp_full.faa",
        rc_core     = _PALM_DIR + "/RdRpCATCH/palm_core.faa",
        lp_core     = _PALM_DIR + "/LucaProt/palm_core.faa",
        rc_extended = _PALM_DIR + "/RdRpCATCH/palm_extended.faa",
        lp_extended = _PALM_DIR + "/LucaProt/palm_extended.faa",
    params:
        outdir         = _PALM_DIR,
        rdrpcatch_glob = _RC_DIR + "/**/*_trimmed_aminoacid_sequences.fasta",
        lucaprot_glob  = _ORF_DIR + "/" + _STEM + "_orfs.faa",
        flank          = 150,
    log:
        out = _ESV_LOG + "/10b_palm_extracted/extract_palm_regions.log",
        err = _ESV_LOG + "/10b_palm_extracted/extract_palm_regions.err",
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
        mkdir -p {params.outdir} $(dirname {log.err})
        python {input.script} \
            --rdrpcatch-motif {input.rdrpcatch_motif} \
            --rdrpcatch-seqs  "{params.rdrpcatch_glob}" \
            --lucaprot-motif  {input.lucaprot_motif} \
            --lucaprot-seqs   "{params.lucaprot_glob}" \
            --outdir          {params.outdir} \
            --flank           {params.flank} \
            > {log.out} 2> {log.err}
        """


# ── Step 10: combine RdRpCATCH + LucaProt into unified sets ──────────────────
# Unlike workflow 03's make_final_output (which also does nr-viral-origin
# filtering via Diamond BLASTp vs nr before combining), esvirtu reference
# sequences are already curated/trusted, so this just concatenates the two
# sources' full/core/extended FASTAs and palm_regions.tsv tables directly --
# no nr-filtering step in this pipeline.

rule combine_rdrp_esvirtu:
    input:
        rc_full     = _PALM_DIR + "/RdRpCATCH/rdrp_full.faa",
        lp_full     = _PALM_DIR + "/LucaProt/rdrp_full.faa",
        rc_core     = _PALM_DIR + "/RdRpCATCH/palm_core.faa",
        lp_core     = _PALM_DIR + "/LucaProt/palm_core.faa",
        rc_extended = _PALM_DIR + "/RdRpCATCH/palm_extended.faa",
        lp_extended = _PALM_DIR + "/LucaProt/palm_extended.faa",
        rc_palm_tsv = _PALM_DIR + "/RdRpCATCH/palm_regions.tsv",
        lp_palm_tsv = _PALM_DIR + "/LucaProt/palm_regions.tsv",
    output:
        combined_full     = _FINAL_DIR + "/combined_full.faa",
        combined_core     = _FINAL_DIR + "/combined_palm_core.faa",
        combined_extended = _FINAL_DIR + "/combined_palm_extended.faa",
        combined_tsv      = _FINAL_DIR + "/combined_palm_regions.tsv",
    log:
        err = _ESV_LOG + "/11_final/combine_rdrp_esvirtu.err",
    threads: config["small_job"]["threads"]
    resources:
        mem_mb_per_cpu  = config["small_job"]["memory"],
        runtime         = config["small_job"]["runtime"],
        cpus_per_task   = config["small_job"]["threads"],
        slurm_partition = config["small_job"]["partition"],
        slurm_account   = config["small_job"]["account"],
    shell:
        """
        mkdir -p {_FINAL_DIR} $(dirname {log.err})
        cat {input.rc_full} {input.lp_full} > {output.combined_full} 2> {log.err}
        cat {input.rc_core} {input.lp_core} > {output.combined_core} 2>> {log.err}
        cat {input.rc_extended} {input.lp_extended} > {output.combined_extended} 2>> {log.err}
        head -n1 {input.rc_palm_tsv} > {output.combined_tsv} 2>> {log.err}
        tail -n +2 {input.rc_palm_tsv} >> {output.combined_tsv} 2>> {log.err}
        tail -n +2 {input.lp_palm_tsv} >> {output.combined_tsv} 2>> {log.err}
        """


# ── Step 11: taxonomy annotation (direct lookup, no DIAMOND needed) ──────────
# ESvirtu reference sequences already have known taxonomy in
# esvirtu_reference.metadata.tsv (per accession) -- no DIAMOND-vs-RVMT search
# is needed like in workflow 03/03c. annotate_esvirtu_taxonomy.py extracts the
# accession from each protein header ("esvirtu_{ACCESSION}") and looks up its
# Phylum/Class/Order/Family/Genus/Species directly. Output TSV matches the
# column shape split_fasta_by_phylum.py expects (qseqid, hit_rank=1, sseqid="",
# Phylum, Class, Order, Family, ...), so that script is reused unchanged.

_ESV_QUERIES = {
    "full_length":   _FINAL_DIR + "/combined_full.faa",
    "palm_core":     _FINAL_DIR + "/combined_palm_core.faa",
    "palm_extended": _FINAL_DIR + "/combined_palm_extended.faa",
}

rule annotate_esvirtu_taxonomy:
    input:
        fasta    = lambda wc: _ESV_QUERIES[wc.query],
        metadata = _ESV_META_TSV,
        script   = "scripts/03_RDRP_identification/annotate_esvirtu_taxonomy.py",
    output:
        tsv = _TAXDIR + "/{query}_taxonomy.tsv",
    wildcard_constraints:
        query = "full_length|palm_core|palm_extended",
    log:
        out = _ESV_LOG + "/12_taxonomy/{query}.log",
        err = _ESV_LOG + "/12_taxonomy/{query}.err",
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
        mkdir -p {_TAXDIR} $(dirname {log.out})
        python {input.script} \
            --fasta    {input.fasta} \
            --metadata {input.metadata} \
            --output   {output.tsv} \
            > {log.out} 2> {log.err}
        """


# ── Step 12: split into per-taxon FASTA files (Phylum / Class / Order / Family)
# full_length is classified directly. palm_core / palm_extended reuse
# full_length's taxonomy assignment for the same seq_id (same header across
# all three protein-region sets, same as workflow 03's Step 22 convention).

rule split_esvirtu_by_phylum:
    input:
        fasta     = _ESV_QUERIES["full_length"],
        annotated = _TAXDIR + "/full_length_taxonomy.tsv",
    output:
        outdir = directory(_PHYLUM_DIR + "/full_length_phylum"),
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
        out = _ESV_LOG + "/13_phylum_cluster/full_length_phylum.log",
        err = _ESV_LOG + "/13_phylum_cluster/full_length_phylum.err",
    shell:
        """
        mkdir -p {output.outdir} $(dirname {log.out})
        python scripts/03_RDRP_identification/split_fasta_by_phylum.py \
            --fasta     {input.fasta} \
            --annotated {input.annotated} \
            --outdir    {output.outdir} \
            > {log.out} 2> {log.err}
        """


rule split_esvirtu_by_phylum_palm:
    input:
        fasta     = lambda wc: _ESV_QUERIES[wc.palmset],
        annotated = _TAXDIR + "/full_length_taxonomy.tsv",
    output:
        outdir = directory(_PHYLUM_DIR + "/{palmset}_phylum"),
    wildcard_constraints:
        palmset = "palm_core|palm_extended",
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
        out = _ESV_LOG + "/13_phylum_cluster/{palmset}_phylum.log",
        err = _ESV_LOG + "/13_phylum_cluster/{palmset}_phylum.err",
    shell:
        """
        mkdir -p {output.outdir} $(dirname {log.out})
        python scripts/03_RDRP_identification/split_fasta_by_phylum.py \
            --fasta     {input.fasta} \
            --annotated {input.annotated} \
            --outdir    {output.outdir} \
            > {log.out} 2> {log.err}
        """


rule split_esvirtu_by_rank:
    input:
        fasta     = _ESV_QUERIES["full_length"],
        annotated = _TAXDIR + "/full_length_taxonomy.tsv",
    output:
        outdir = directory(_PHYLUM_DIR + "/full_length_{rank}"),
    wildcard_constraints:
        rank = "class|order|family",
    params:
        rank_col = lambda wc: wc.rank.capitalize(),
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
        out = _ESV_LOG + "/13_phylum_cluster/full_length_{rank}.log",
        err = _ESV_LOG + "/13_phylum_cluster/full_length_{rank}.err",
    shell:
        """
        mkdir -p {output.outdir} $(dirname {log.out})
        python scripts/03_RDRP_identification/split_fasta_by_phylum.py \
            --fasta     {input.fasta} \
            --annotated {input.annotated} \
            --outdir    {output.outdir} \
            --rank      {params.rank_col} \
            > {log.out} 2> {log.err}
        """


rule split_esvirtu_by_rank_palm:
    input:
        fasta     = lambda wc: _ESV_QUERIES[wc.palmset],
        annotated = _TAXDIR + "/full_length_taxonomy.tsv",
    output:
        outdir = directory(_PHYLUM_DIR + "/{palmset}_{rank}"),
    wildcard_constraints:
        palmset = "palm_core|palm_extended",
        rank    = "class|order|family",
    params:
        rank_col = lambda wc: wc.rank.capitalize(),
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
        out = _ESV_LOG + "/13_phylum_cluster/{palmset}_{rank}.log",
        err = _ESV_LOG + "/13_phylum_cluster/{palmset}_{rank}.err",
    shell:
        """
        mkdir -p {output.outdir} $(dirname {log.out})
        python scripts/03_RDRP_identification/split_fasta_by_phylum.py \
            --fasta     {input.fasta} \
            --annotated {input.annotated} \
            --outdir    {output.outdir} \
            --rank      {params.rank_col} \
            > {log.out} 2> {log.err}
        """
