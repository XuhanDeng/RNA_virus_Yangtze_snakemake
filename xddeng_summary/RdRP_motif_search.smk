"""
RdRP motif A-D identification pipeline (Snakemake)
Based on RVMT (Neri et al.) — all commands inline, no wrapper scripts needed.

Usage:
    snakemake -s RdRP_motif_search.smk --configfile config.yaml -j 16

Required config keys (see config.yaml):
    input_faa        : path to your candidate RdRP protein sequences (.faa)
    sequence_library : path to Zenodo Sequence_Library/ directory
                       (must contain mot.1/ mot.2/ mot.3/ mot.4/ subdirs)
    outdir           : base output directory
    threads          : number of CPU threads
    prefix           : label for output files (e.g. "MyRun")
"""

configfile: "config.yaml"

import os

MOTIFS  = [1, 2, 3, 4]
OUTDIR  = config["outdir"]
FAA     = config["input_faa"]
SEQ_LIB = config["sequence_library"]
THREADS = config["threads"]
PREFIX  = config["prefix"]
STEM    = os.path.basename(FAA).replace(".faa", "")


# ─────────────────────────────────────────────────────────────────────────────
# Target rule
# ─────────────────────────────────────────────────────────────────────────────
rule all:
    input:
        # Profiles
        expand("{outdir}/profiles/mot.{motif}/HMMfiles/profiles.hmm",
               outdir=OUTDIR, motif=MOTIFS),
        expand("{outdir}/profiles/mot.{motif}/HHMfiles/db",
               outdir=OUTDIR, motif=MOTIFS),
        expand("{outdir}/profiles/mot.{motif}/MMseqs2_profiles/MM",
               outdir=OUTDIR, motif=MOTIFS),
        # Search results
        expand("{outdir}/search/{tool}/mot.{motif}.tsv",
               outdir=OUTDIR,
               tool=["hhalign", "hmmsearch", "psiblast", "mmseqs", "diamond"],
               motif=MOTIFS),


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1a  Add consensus sequence to each MSA
#          hhconsensus -M 50 -cov 50
#          Input : Sequence_Library/mot.{motif}/*_woCon.afa
#          Output: outdir/profiles/mot.{motif}/msaFiles/*.Cons.msa.afa
# ─────────────────────────────────────────────────────────────────────────────
rule add_consensus:
    input:
        msa_dir = SEQ_LIB + "/mot.{motif}"
    output:
        done = "{outdir}/profiles/mot.{motif}/msaFiles/.consensus_done"
    params:
        out_msa = "{outdir}/profiles/mot.{motif}/msaFiles"
    threads: THREADS
    log:
        "{outdir}/logs/consensus_mot{motif}.log"
    shell:
        """
        mkdir -p {params.out_msa}

        # Copy original MSAs into msaFiles/
        cp {input.msa_dir}/*.afa {params.out_msa}/

        # Add consensus sequence on top of each alignment
        # -M 50: use columns with >50% occupancy for consensus
        # -cov 50: minimum 50% coverage threshold
        parallel -j {threads} \
            hhconsensus -M 50 -cov 50 \
                -i {{}} \
                -oa2m {params.out_msa}/{{/.}}.Cons.msa.afa \
            ::: {params.out_msa}/*.afa \
        >> {log} 2>&1

        # Remove the first line (hhconsensus inserts a blank/comment line)
        parallel -j {threads} sed -i '1d' {{}} \
            ::: {params.out_msa}/*.Cons.msa.afa \
        >> {log} 2>&1

        touch {output.done}
        """


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1b  Build HMMER3 profile database
#          hmmbuild per cluster → cat into one .hmm file
#          Input : *.Cons.msa.afa (from add_consensus)
#          Output: outdir/profiles/mot.{motif}/HMMfiles/profiles.hmm
# ─────────────────────────────────────────────────────────────────────────────
rule build_hmm:
    input:
        done = "{outdir}/profiles/mot.{motif}/msaFiles/.consensus_done"
    output:
        hmm_db = "{outdir}/profiles/mot.{motif}/HMMfiles/profiles.hmm"
    params:
        msa_dir  = "{outdir}/profiles/mot.{motif}/msaFiles",
        hmm_dir  = "{outdir}/profiles/mot.{motif}/HMMfiles",
        log_dir  = "{outdir}/profiles/mot.{motif}/logs"
    threads: THREADS
    log:
        "{outdir}/logs/hmmbuild_mot{motif}.log"
    shell:
        """
        mkdir -p {params.hmm_dir} {params.log_dir}

        # Build one HMM profile per consensus MSA, in parallel
        # --informat a2m : input is a2m/afa format
        # -n {/.}        : profile name = filename without extension
        parallel -j {threads} \
            hmmbuild --informat a2m \
                -n {{/.}} \
                -o {params.log_dir}/{{/.}}.hmm.log \
                {params.hmm_dir}/{{/.}}.hmm \
                {{}} \
            ::: {params.msa_dir}/*.Cons.msa.afa \
        >> {log} 2>&1

        # Concatenate all per-profile .hmm files into one searchable database
        cat {params.hmm_dir}/*.hmm > {output.hmm_db}
        """


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1c  Build HH-suite profile database (ffindex)
#          hhmake per cluster → ffindex_build
#          Input : *.Cons.msa.afa
#          Output: outdir/profiles/mot.{motif}/HHMfiles/db  (ffindex prefix)
# ─────────────────────────────────────────────────────────────────────────────
rule build_hhm:
    input:
        done = "{outdir}/profiles/mot.{motif}/msaFiles/.consensus_done"
    output:
        db = "{outdir}/profiles/mot.{motif}/HHMfiles/db"
    params:
        msa_dir = "{outdir}/profiles/mot.{motif}/msaFiles",
        hhm_dir = "{outdir}/profiles/mot.{motif}/HHMfiles"
    threads: THREADS
    log:
        "{outdir}/logs/hhmake_mot{motif}.log"
    shell:
        """
        mkdir -p {params.hhm_dir}

        # Build one .hhm profile per consensus MSA
        # -M a2m : treat input as a2m format
        # -v 2   : verbosity
        parallel -j {threads} \
            hhmake -v 2 \
                -name {{/.}} \
                -i {{}} \
                -o {params.hhm_dir}/{{/.}}.hhm \
                -M a2m \
            ::: {params.msa_dir}/*.Cons.msa.afa \
        >> {log} 2>&1

        # Index all .hhm files into a single ffindex database
        # This produces db (data file) and db.index
        ffindex_build \
            {params.hhm_dir}/db \
            {params.hhm_dir}/db.index \
            {params.hhm_dir}/*.hhm \
        >> {log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1d  Build MMseqs2 profile database
#          Converts the ffindex HHM database into MMseqs2 profile format
#          Input : HHMfiles/db  (from build_hhm)
#          Output: outdir/profiles/mot.{motif}/MMseqs2_profiles/MM
# ─────────────────────────────────────────────────────────────────────────────
rule build_mmseqs_profile:
    input:
        hhm_db = "{outdir}/profiles/mot.{motif}/HHMfiles/db"
    output:
        mm_db = "{outdir}/profiles/mot.{motif}/MMseqs2_profiles/MM"
    params:
        mm_dir  = "{outdir}/profiles/mot.{motif}/MMseqs2_profiles",
        hhm_db  = "{outdir}/profiles/mot.{motif}/HHMfiles/db"
    log:
        "{outdir}/logs/mmseqs_profile_mot{motif}.log"
    shell:
        """
        mkdir -p {params.mm_dir}

        # Convert HH-suite ffindex profile DB into MMseqs2 profile DB
        mmseqs convertprofiledb \
            {params.hhm_db} \
            {output.mm_db} \
        >> {log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# STEP 1e  Collect singleton sequences
#          MSAs with only 1 sequence — used as Diamond BLASTp query
#          Input : msaFiles/*.afa
#          Output: outdir/profiles/mot.{motif}/singletons.faa
# ─────────────────────────────────────────────────────────────────────────────
rule collect_singletons:
    input:
        done = "{outdir}/profiles/mot.{motif}/msaFiles/.consensus_done"
    output:
        singletons = "{outdir}/profiles/mot.{motif}/singletons.faa"
    params:
        msa_dir = "{outdir}/profiles/mot.{motif}/msaFiles"
    log:
        "{outdir}/logs/singletons_mot{motif}.log"
    shell:
        """
        # Find all original .afa files (not Cons) with exactly 1 sequence
        # and concatenate them as the singleton sequence set
        grep -l "" {params.msa_dir}/*.afa \
            | while read f; do
                n=$(grep -c "^>" "$f" || true)
                if [ "$n" -eq 1 ]; then
                    cat "$f"
                fi
              done \
            > {output.singletons} \
        2>> {log}

        # Ensure file exists even if no singletons found
        touch {output.singletons}
        """


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2a  HHsearch  (profile–profile, highest sensitivity)
#          Query  : split your RdRPs into per-sequence .faa, run hhsearch each
#          Subject: HH-suite profile DB (HHMfiles/db)
#          Output : outdir/search/hhalign/mot.{motif}.tsv
#
#          Output columns (13, HHsearch blasttab):
#            1  subject_name   matched motif profile
#            2  pCoverage      (p2-p1+1)/pL
#            3  ali_len        alignment length
#            4  pL             profile length
#            5  mismatch
#            6  gapOpen
#            7  q1             start on YOUR sequence (AA)
#            8  q2             end on YOUR sequence (AA)
#            9  p1             start on motif profile
#           10  p2             end on motif profile
#           11  Probab         HHsearch probability
#           12  evalue
#           13  score
# ─────────────────────────────────────────────────────────────────────────────
rule search_hhalign:
    input:
        faa   = FAA,
        hhmdb = "{outdir}/profiles/mot.{motif}/HHMfiles/db"
    output:
        tsv = "{outdir}/search/hhalign/mot.{motif}.tsv"
    params:
        split_dir  = "{outdir}/search/hhalign/mot.{motif}_split",
        result_dir = "{outdir}/search/hhalign/mot.{motif}_results"
    threads: THREADS
    log:
        "{outdir}/logs/hhalign_mot{motif}.log"
    shell:
        """
        mkdir -p {params.split_dir} {params.result_dir}

        # Split multi-FASTA into one file per sequence
        # splitfasta.pl is part of HH-suite
        cd {params.split_dir}
        splitfasta.pl {input.faa} -ext .faa >> {log} 2>&1
        cd -

        # Run hhsearch on each sequence in parallel
        # -M 50     : minimum column occupancy 50% for consensus
        # -e 0.5    : E-value cutoff (loose; R parsing applies stricter filter)
        # -hide_cons: suppress consensus row in output
        # -Z 100000 : return up to 100,000 hits
        ls {params.split_dir}/*.faa \
            | parallel -j {threads} \
                hhsearch \
                    -M 50 \
                    -e 0.5 \
                    -i {{}} \
                    -d {input.hhmdb} \
                    -o /dev/null \
                    -cpu 1 \
                    -hide_cons \
                    -Z 100000 \
                    -blasttab {params.result_dir}/{{/.}}_out.tsv \
            >> {log} 2>&1

        # Concatenate all per-sequence results
        cat {params.result_dir}/*.tsv > {output.tsv} 2>> {log}

        # Remove consensus suffix that hhsearch inserts into profile names
        sed -i 's|.Cons.msa||g' {output.tsv}
        """


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2b  HMMER hmmsearch
#          Query  : your RdRPs (.faa)
#          Subject: HMMER3 profile DB (HMMfiles/profiles.hmm)
#          Output : outdir/search/hmmsearch/mot.{motif}.tsv
#
#          Output columns (23, HMMER domtblout, comment lines stripped):
#            1  target_name    your sequence ID
#            2  target_acc     "-"
#            3  qL             your sequence length
#            4  query_name     motif profile name
#            5  query_acc      "-"
#            6  pL             motif profile length
#            7  E-value        full-sequence E-value
#            8  score          full-sequence bit score
#            9  bias
#           10  #              domain number
#           11  of             total domains
#           12  c-Evalue       conditional E-value
#           13  i-Evalue       independent E-value
#           14  score2         domain bit score
#           15  bias2
#           16  p1             start on motif profile
#           17  p2             end on motif profile
#           18  q1             start on YOUR sequence (AA)
#           19  q2             end on YOUR sequence (AA)
#           20  env_from
#           21  env_to
#           22  acc            mean posterior probability
#           23  description
# ─────────────────────────────────────────────────────────────────────────────
rule search_hmmsearch:
    input:
        faa = FAA,
        hmm = "{outdir}/profiles/mot.{motif}/HMMfiles/profiles.hmm"
    output:
        tsv = "{outdir}/search/hmmsearch/mot.{motif}.tsv"
    params:
        raw = "{outdir}/search/hmmsearch/mot.{motif}_raw.tsv"
    threads: THREADS
    log:
        "{outdir}/logs/hmmsearch_mot{motif}.log"
    shell:
        """
        mkdir -p {OUTDIR}/search/hmmsearch

        # --noali  : suppress alignment output (faster)
        # -E 0.5   : E-value cutoff (loose; R parsing applies stricter filter)
        # --domtblout: per-domain tabular output (gives q1/q2/p1/p2 per hit)
        hmmsearch \
            --noali \
            --cpu {threads} \
            -E 0.5 \
            --incE 0.5 \
            --domtblout {params.raw} \
            {input.hmm} \
            {input.faa} \
        >> {log} 2>&1

        # Strip comment lines starting with #
        grep -v "^#" {params.raw} > {output.tsv} 2>> {log} || true
        """


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2c  PSI-BLAST
#          Query  : motif consensus MSAs (*Cons.msa.afa in msaFiles/)
#          Subject: your RdRPs (built as BLAST DB)
#          Output : outdir/search/psiblast/mot.{motif}.tsv
#
#          Output columns (12):
#            1  sseqid       your sequence ID  (= RdRp_id)
#            2  pident       % identity
#            3  sstart       start on YOUR sequence  (= q1)
#            4  send         end on YOUR sequence    (= q2)
#            5  qstart       start on motif MSA      (= p1)
#            6  qend         end on motif MSA        (= p2)
#            7  slen         your sequence length    (= qL)
#            8  qlen         motif MSA length        (= pL)
#            9  length       alignment length        (= ali_len)
#           10  evalue
#           11  bitscore
#           12  profile_name  motif cluster name (appended from filename)
#
#          NOTE: query/subject are REVERSED vs other tools.
#                Motif MSA = query; your RdRP = subject.
#                So q1/q2 are in sstart/send (cols 3/4), not qstart/qend.
# ─────────────────────────────────────────────────────────────────────────────
rule search_psiblast:
    input:
        faa       = FAA,
        consensus = "{outdir}/profiles/mot.{motif}/msaFiles/.consensus_done"
    output:
        tsv = "{outdir}/search/psiblast/mot.{motif}.tsv"
    params:
        db_dir  = "{outdir}/search/psiblast/mot.{motif}_blastdb",
        tmp_dir = "{outdir}/search/psiblast/mot.{motif}_tmp",
        msa_dir = "{outdir}/profiles/mot.{motif}/msaFiles"
    threads: THREADS
    log:
        "{outdir}/logs/psiblast_mot{motif}.log"
    shell:
        """
        mkdir -p {params.db_dir} {params.tmp_dir}

        # Build BLAST protein database from your RdRP sequences
        makeblastdb \
            -in {input.faa} \
            -dbtype prot \
            -title {STEM} \
            -out {params.db_dir}/{STEM} \
        >> {log} 2>&1

        # Run PSI-BLAST for each consensus MSA cluster against your RdRP DB
        # -word_size 2       : shorter words for more sensitivity
        # -threshold 9       : neighbourhood word score threshold
        # -dbsize 20000000   : normalise E-values to fixed DB size
        # -ignore_msa_master : don't use first row as query sequence
        # -qcov_hsp_perc 1   : require at least 1% query coverage
        for cluster in {params.msa_dir}/*.Cons.msa.afa; do
            cluster_name=$(basename "$cluster" ".msa.Cons.msa.afa")
            psiblast \
                -word_size 2 \
                -evalue 0.5 \
                -max_target_seqs 100000 \
                -threshold 9 \
                -dbsize 20000000 \
                -ignore_msa_master \
                -in_msa "$cluster" \
                -qcov_hsp_perc 1 \
                -num_threads {threads} \
                -db {params.db_dir}/{STEM} \
                -outfmt "6 sseqid pident sstart send qstart qend slen qlen length evalue bitscore" \
                -out {params.tmp_dir}/"${{cluster_name}}"_out.tsv \
            >> {log} 2>&1

            # Append cluster name as final column
            awk -v cls="$cluster_name" \
                'BEGIN{{OFS="\t"}} {{print $0, cls}}' \
                {params.tmp_dir}/"${{cluster_name}}"_out.tsv \
            >> {output.tsv}
        done
        """


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2d  MMseqs2 profile search
#          Query  : your RdRPs (converted to MMseqs2 sequence DB)
#          Subject: MMseqs2 profile DB (MMseqs2_profiles/MM)
#          Output : outdir/search/mmseqs/mot.{motif}.tsv
#
#          Output columns (19):
#            1  query       your sequence ID
#            2  target      motif profile cluster name
#            3  evalue
#            4  gapopen
#            5  pident      % identity
#            6  nident      number of identical residues
#            7  qstart      start on YOUR sequence  (= q1)
#            8  qend        end on YOUR sequence    (= q2)
#            9  qlen        your sequence length    (= qL)
#           10  tstart      start on motif profile  (= p1)
#           11  tend        end on motif profile    (= p2)
#           12  tlen        motif profile length    (= pL)
#           13  alnlen      alignment length
#           14  raw         raw score
#           15  bits        bit score
#           16  qframe      reading frame (0 for protein)
#           17  mismatch
#           18  qcov        query coverage
#           19  tcov        target (profile) coverage
# ─────────────────────────────────────────────────────────────────────────────
rule search_mmseqs:
    input:
        faa  = FAA,
        mmdb = "{outdir}/profiles/mot.{motif}/MMseqs2_profiles/MM"
    output:
        tsv = "{outdir}/search/mmseqs/mot.{motif}.tsv"
    params:
        query_db   = "{outdir}/search/mmseqs/mot.{motif}_querydb/DB",
        result_db  = "{outdir}/search/mmseqs/mot.{motif}_results/DB",
        tmp_dir    = "{outdir}/search/mmseqs/mot.{motif}_tmp"
    threads: THREADS
    log:
        "{outdir}/logs/mmseqs_mot{motif}.log"
    shell:
        """
        mkdir -p \
            $(dirname {params.query_db}) \
            $(dirname {params.result_db}) \
            {params.tmp_dir}

        # Convert your .faa to MMseqs2 sequence database
        mmseqs createdb {input.faa} {params.query_db} >> {log} 2>&1

        # Profile search
        # -k 6    : k-mer size 6 (shorter = more sensitive)
        # -s 7.5  : sensitivity 7.5 (max sensitivity)
        # -e 0.5  : E-value cutoff (loose)
        # --max-seqs 6000 : return up to 6000 hits per query
        mmseqs search \
            {params.query_db} \
            {input.mmdb} \
            {params.result_db} \
            {params.tmp_dir} \
            -k 6 \
            -s 7.5 \
            -e 0.5 \
            --threads {threads} \
            --split-memory-limit 100G \
            --max-seqs 6000 \
        >> {log} 2>&1

        # Convert result DB to TSV
        mmseqs convertalis \
            {params.query_db} \
            {input.mmdb} \
            {params.result_db} \
            {output.tsv} \
            --format-output \
                "query,target,evalue,gapopen,pident,nident,qstart,qend,qlen,tstart,tend,tlen,alnlen,raw,bits,qframe,mismatch,qcov,tcov" \
        >> {log} 2>&1
        """


# ─────────────────────────────────────────────────────────────────────────────
# STEP 2e  Diamond BLASTp  (strictest at search time: id≥60%, qcov≥50%)
#          Query  : motif singleton sequences (singletons.faa)
#          Subject: your RdRPs (built as Diamond DB)
#          Output : outdir/search/diamond/mot.{motif}.tsv
#
#          Output columns (14):
#            1  sseqid    your sequence ID         (= RdRp_id)
#            2  qseqid    motif singleton ID, prefixed "mot.{motif}."
#            3  evalue
#            4  gapopen
#            5  pident    % identity
#            6  length    alignment length         (= ali_len)
#            7  qstart    start on singleton       (= p1)
#            8  qend      end on singleton         (= p2)
#            9  sstart    start on YOUR sequence   (= q1)
#           10  send      end on YOUR sequence     (= q2)
#           11  qlen      singleton length         (= pL)
#           12  slen      your sequence length     (= qL)
#           13  mismatch
#           14  bitscore
#
#          NOTE: query/subject are REVERSED.
#                Singletons = query; your RdRPs = subject.
#                q1/q2 are in cols 9/10 (sstart/send), NOT 7/8.
# ─────────────────────────────────────────────────────────────────────────────
rule search_diamond:
    input:
        faa        = FAA,
        singletons = "{outdir}/profiles/mot.{motif}/singletons.faa"
    output:
        tsv = "{outdir}/search/diamond/mot.{motif}.tsv"
    params:
        db_dir = "{outdir}/search/diamond/mot.{motif}_db",
        raw    = "{outdir}/search/diamond/mot.{motif}_raw.tsv"
    threads: THREADS
    log:
        "{outdir}/logs/diamond_mot{motif}.log"
    shell:
        """
        mkdir -p {params.db_dir}

        # Build Diamond database from your RdRP sequences (subject)
        diamond makedb \
            --in {input.faa} \
            -d {params.db_dir}/{STEM} \
            -p {threads} \
        >> {log} 2>&1

        # BLASTp: singleton sequences (query) vs your RdRPs (subject)
        # --query-cover 50 : require 50% of the singleton to align
        # --id 60          : require 60% identity
        # --more-sensitive : slower but catches more divergent hits
        # -k 0             : return all hits (no limit)
        # -b 3.0           : larger memory block for speed
        diamond blastp \
            -q {input.singletons} \
            --db {params.db_dir}/{STEM} \
            -p {threads} \
            -b 3.0 \
            -k 0 \
            --query-cover 50 \
            --id 60 \
            --more-sensitive \
            --evalue 0.05 \
            --outfmt 6 sseqid qseqid evalue gapopen pident length \
                        qstart qend sstart send qlen slen mismatch bitscore \
            -o {params.raw} \
        >> {log} 2>&1

        # Prepend "mot.{motif}." to the query (singleton) ID column
        awk -F'\t' -vOFS='\t' -v m={wildcards.motif} \
            '{{ $2 = "mot." m "." $2 }}1' \
            {params.raw} > {output.tsv}
        """
