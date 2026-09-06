# ============================================================
# 00_lucaprot.smk
# Run on login node — requires internet and conda
# ============================================================

# dry-run
snakemake --snakefile workflow/00_lucaprot.smk --cores 1 -n -p

snakemake --snakefile workflow/00_lucaprot.smk --cores 1  --use-conda

# unlock
snakemake --snakefile workflow/00_lucaprot.smk --unlock
nohup snakemake --snakefile workflow/00_lucaprot.smk --sdm conda --conda-create-envs-only \
    > log/00_lucaprot/snakemake_envs.log 2>&1 &
echo "PID: $!"

# formal run (login node, no SLURM needed)
mkdir -p log/00_lucaprot
nohup snakemake --snakefile workflow/00_lucaprot.smk \
    --cores 1 \
    --rerun-triggers input \
    --use-conda \
    > log/00_lucaprot/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 00_setup_RNA_database.smk
# Run on login node — HPC compute nodes lack internet access
# ============================================================

# dry-run
snakemake --snakefile workflow/00_setup_RNA_database.smk --cores 2 --use-conda -n -p --rerun-triggers input

# build conda envs only
snakemake --snakefile workflow/00_setup_RNA_database.smk --sdm conda --conda-create-envs-only



nohup snakemake --snakefile workflow/00_setup_RNA_database.smk --cores 2 --use-conda --rerun-triggers input > log/00_setup_RNA_database/snakemake.log 2>&1 &




# unlock
snakemake --snakefile workflow/00_setup_RNA_database.smk --unlock

# formal run
snakemake --snakefile workflow/00_setup_RNA_database.smk --cores 2 --use-conda --rerun-triggers input


# ============================================================
# 01_bacteria_genome_assembly.smk
# ============================================================

# dry-run
snakemake --snakefile workflow/01_bacteria_genome_assembly.smk --use-conda --cores 1 --rerun-triggers input -n -p 

# build conda envs only
snakemake --snakefile workflow/01_bacteria_genome_assembly.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/01_bacteria_genome_assembly.smk --unlock

# formal run
mkdir -p log/01_bacteria_genome_assembly
nohup snakemake --snakefile workflow/01_bacteria_genome_assembly.smk \
    --executor slurm \
    --jobs 3 \
    --use-conda \
    --retries 1 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/01_bacteria_genome_assembly/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 02_RNA_virus_assembly.smk
# ============================================================

# dry-run
snakemake --snakefile workflow/02_RNA_virus_assembly.smk --use-conda --cores 1 -n -p

# build conda envs only
snakemake --snakefile workflow/02_RNA_virus_assembly.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/02_RNA_virus_assembly.smk --unlock

# formal run
mkdir -p log/02_RNA_virus_assembly
nohup snakemake --snakefile workflow/02_RNA_virus_assembly.smk \
    --executor slurm \
    --jobs 50 \
    --use-conda \
    --retries 1 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/02_RNA_virus_assembly/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 02_RNA_virus_assembly_sortmerna.smk
# ============================================================

# dry-run
snakemake --snakefile workflow/02_RNA_virus_assembly_sortmerna.smk --use-conda --cores 1 -n -p --rerun-triggers input  --until sortmerna_rrna_removal

# build conda envs only
snakemake --snakefile workflow/02_RNA_virus_assembly_sortmerna.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/02_RNA_virus_assembly_sortmerna.smk --unlock

# formal run
mkdir -p log/02_RNA_virus_assembly
nohup snakemake --snakefile workflow/02_RNA_virus_assembly_sortmerna.smk \
    --executor slurm \
    --jobs 16 \
    --use-conda \
    --retries 1 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/02_RNA_virus_assembly/snakemake_sortmerna.log 2>&1 &
echo "PID: $!"

mkdir -p log/02_RNA_virus_assembly
nohup snakemake --snakefile workflow/02_RNA_virus_assembly_sortmerna.smk \
    --executor slurm --jobs 50 --use-conda --retries 1 \
    --printshellcmds --slurm-no-account --rerun-triggers input \
    --latency-wait 60 \
    --until sortmerna_rrna_removal \
    > log/02_RNA_virus_assembly/snakemake_sortmerna.log 2>&1 &
echo "PID: $!"


# ============================================================
# 03_RDRP_identification.smk
# ============================================================


# dry-run (ICTV enabled)
snakemake --snakefile workflow/03_RDRP_identification.smk --use-conda --cores 4     --rerun-triggers input -n -p

snakemake --snakefile workflow/03_RDRP_identification.smk --use-conda --cores 4  --rerun-triggers input -n -p --until motif_search_diamond



# build conda envs only
snakemake --snakefile workflow/03_RDRP_identification.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/03_RDRP_identification.smk --unlock

# formal run (ICTV disabled)
mkdir -p log/03_RDRP_identification
nohup snakemake --snakefile workflow/03_RDRP_identification.smk \
    --executor slurm \
    --jobs 64 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/03_RDRP_identification/snakemake.log 2>&1 &
echo "PID: $!"

# formal run (ICTV enabled)
mkdir -p log/03_RDRP_identification
    nohup snakemake --snakefile workflow/03_RDRP_identification.smk \
    --executor slurm \
    --jobs 64 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/03_RDRP_identification/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 03c_ref_palm_annot.smk
# Run after 03_RDRP_identification.smk (through Step 22, split_by_phylum /
# split_by_rank) is complete -- {taxon}_ref.faa files must already exist for
# all 4 ranks (phylum/class/order/family).
# ============================================================

# dry-run
snakemake --snakefile workflow/03c_ref_palm_annot.smk --use-conda --cores 4 --rerun-triggers input -n -p

# unlock
snakemake --snakefile workflow/03c_ref_palm_annot.smk --unlock

# formal run
mkdir -p log/03c_ref_palm_annot
nohup snakemake --snakefile workflow/03c_ref_palm_annot.smk \
    --executor slurm \
    --jobs 60 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    --keep-going \
    > log/03c_ref_palm_annot/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 04_modified_esvirtue.smk
# ============================================================

# dry-run
snakemake --snakefile workflow/04_modified_esvirtue.smk --use-conda --cores 1 -n -p

# build conda envs only
snakemake --snakefile workflow/04_modified_esvirtue.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/04_modified_esvirtue.smk --unlock

# formal run
mkdir -p log/04_modified_esvirtue
nohup snakemake --snakefile workflow/04_modified_esvirtue.smk \
    --executor slurm \
    --jobs 50 \
    --use-conda \
    --retries 1 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/04_modified_esvirtue/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 04_modified_esvirtue_ribodetector.smk
# ============================================================

# dry-run
snakemake --snakefile workflow/04_modified_esvirtue_ribodetector.smk --use-conda --cores 1 -n -p   --rerun-triggers input
#run in login
snakemake --snakefile workflow/04_modified_esvirtue_ribodetector.smk --use-conda --cores 2 --rerun-triggers input
# build conda envs only
snakemake --snakefile workflow/04_modified_esvirtue_ribodetector.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/04_modified_esvirtue_ribodetector.smk --unlock

# formal run
mkdir -p log/04_modified_esvirtue_ribodetector
nohup snakemake --snakefile workflow/04_modified_esvirtue_ribodetector.smk \
    --executor slurm \
    --jobs 50 \
    --use-conda \
    --retries 1 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/04_modified_esvirtue_ribodetector/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 03b_esvirtu_rdrp_identification.smk
# Run after 04_modified_esvirtue_ribodetector.smk is complete
# ============================================================

# dry-run
snakemake --snakefile workflow/03b_esvirtu_rdrp_identification.smk --use-conda --cores 1 --rerun-triggers input -n -p


snakemake --snakefile workflow/03b_esvirtu_rdrp_identification.smk --use-conda --cores 4 --rerun-triggers input --until extract_esvirtu_sequences


# build conda envs only
snakemake --snakefile workflow/03b_esvirtu_rdrp_identification.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/03b_esvirtu_rdrp_identification.smk --unlock

# formal run
mkdir -p log/03b_esvirtu_rdrp_identification
nohup snakemake --snakefile workflow/03b_esvirtu_rdrp_identification.smk \
    --executor slurm \
    --jobs 16 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/03b_esvirtu_rdrp_identification/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 05_esvirtue_correlation_analysis.smk
# Run after 04_modified_esvirtue*.smk workflows are complete
# ============================================================

# dry-run
snakemake --snakefile workflow/05_esvirtu_correlation_analysis.smk --use-conda --cores 2 --rerun-triggers input -n -p 

# build conda envs only
snakemake --snakefile workflow/05_esvirtu_correlation_analysis.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/05_esvirtu_correlation_analysis.smk --unlock

# formal run
mkdir -p log/05_esvirtu_correlation_analysis
nohup snakemake --snakefile workflow/05_esvirtu_correlation_analysis.smk \
    --executor slurm \
    --jobs 8 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/05_esvirtu_correlation_analysis/snakemake.log 2>&1 &
echo "PID: $!"

# ============================================================
# 101_kraken_check.smk
# ============================================================

# dry-run
snakemake --snakefile workflow/101_kraken_check.smk --use-conda --cores 1 -n -p

# build conda envs only
snakemake --snakefile workflow/101_kraken_check.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/101_kraken_check.smk --unlock

# formal run
mkdir -p log/101_kraken_check
nohup snakemake --snakefile workflow/101_kraken_check.smk \
    --executor slurm \
    --jobs 50 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/101_kraken_check/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 07_RDRP_tree_construction.smk
# Run after 03_RDRP_identification.smk is complete (Steps 1–8 + ICTV)
# ============================================================

# dry-run
snakemake --snakefile workflow/07_RDRP_tree_construction.smk --use-conda --cores 2 --rerun-triggers input -n -p 

# build conda envs only
snakemake --snakefile workflow/07_RDRP_tree_construction.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/07_RDRP_tree_construction.smk --unlock

# formal run
mkdir -p log/07_RDRP_tree_construction
nohup snakemake --snakefile workflow/07_RDRP_tree_construction.smk \
    --executor slurm \
    --jobs 10 \
    --use-conda \
    --retries 2 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/07_RDRP_tree_construction/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 100_figure.smk
# ============================================================

# dry-run
snakemake --snakefile workflow/100_figure.smk --use-conda --cores 2 --rerun-triggers input -n -p

# build conda envs only
snakemake --snakefile workflow/100_figure.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/100_figure.smk --unlock

# formal run
mkdir -p log/100_network_figure/1_spearman_network_figure
nohup snakemake --snakefile workflow/100_figure.smk \
    --executor slurm \
    --jobs 32 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/100_network_figure/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 06_contig_bait.smk
# Run after 03_RDRP_identification.smk is complete
# ============================================================

# dry-run
snakemake --snakefile workflow/06_contig_bait.smk --use-conda --cores 1 --rerun-triggers input -n -p

# build conda envs only
snakemake --snakefile workflow/06_contig_bait.smk --sdm conda --conda-create-envs-only

# dry-run
snakemake --snakefile workflow/06_contig_bait.smk --use-conda --cores 1 --rerun-triggers input --until extract_rdrp_bait_contigs

# unlock
snakemake --snakefile workflow/06_contig_bait.smk --unlock

# formal run
mkdir -p log/06_contig_bait
nohup snakemake --snakefile workflow/06_contig_bait.smk \
    --executor slurm \
    --jobs 8 \
    --use-conda \
    --retries 1 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/06_contig_bait/snakemake.log 2>&1 &
echo "PID: $!"





snakemake --snakefile workflow/06_contig_bait_memory.smk --use-conda --cores 4 --rerun-triggers input -n -p

snakemake --snakefile workflow/06_contig_bait_memory.smk --unlock


nohup snakemake --snakefile workflow/06_contig_bait_memory.smk \
    --executor slurm \
    --jobs 8 \
    --use-conda \
    --retries 1 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/06_contig_bait/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 08_ssDNA_identification.smk
# Run after 02_RNA_virus_assembly.smk is complete
# ============================================================

# dry-run
snakemake --snakefile workflow/08_other_virus_identification.smk --use-conda --cores 4 --rerun-triggers input -n -p

# build conda envs only
snakemake --snakefile workflow/08_other_virus_identification.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/08_other_virus_identification.smk --unlock

# formal run
mkdir -p log/08_other_virus_identification
nohup snakemake --snakefile workflow/08_other_virus_identification.smk \
    --executor slurm \
    --jobs 50 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 15 \
    --rerun-incomplete \
    > log/08_other_virus_identification/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 99_tables.smk
# Run on login node after all upstream workflows are complete
# ============================================================

# dry-run
snakemake --snakefile workflow/99_tables.smk --use-conda --cores 3 --rerun-triggers input -n -p 

snakemake --snakefile workflow/99_tables.smk --use-conda --cores 3 --rerun-triggers  input -n -p 


snakemake --snakefile workflow/99_tables.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/99_tables.smk --unlock

# formal run
mkdir -p log/99_tables
nohup snakemake --snakefile workflow/99_tables.smk \
    --executor slurm \
    --jobs 50 \
    --use-conda \
    --retries 1 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 15 \
    --rerun-incomplete \
    --keep-going\
    > log/99_tables/snakemake.log 2>&1 &
echo "PID: $!"
echo "PID: $!"


# ============================================================
# 101_cluster_RDRP_tree_construction.smk
# Run after 03_RDRP_identification.smk and 99_tables.smk are complete
# ============================================================

# dry-run
snakemake --snakefile workflow/101_cluster_RDRP_tree_construction.smk --use-conda --cores 2 --rerun-triggers input -n -p


snakemake --snakefile workflow/101_cluster_RDRP_tree_construction.smk --use-conda --cores 2 --rerun-triggers input --until itol_colorstrip



# build conda envs only
snakemake --snakefile workflow/101_cluster_RDRP_tree_construction.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/101_cluster_RDRP_tree_construction.smk --unlock

# formal run
mkdir -p log/101_cluster_RDRP_tree_construction
mkdir -p result/101_cluster_RDRP_tree_construction
nohup snakemake --snakefile workflow/101_cluster_RDRP_tree_construction.smk \
    --executor slurm \
    --jobs 10 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    --keep-going\
    > log/101_cluster_RDRP_tree_construction/snakemake.log 2>&1 &
echo "PID: $!"


nohup snakemake --snakefile workflow/101_cluster_RDRP_tree_construction.smk \
    --executor slurm \
    --jobs 10 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    --until rename_tree_tips_ictv\
    > log/101_cluster_RDRP_tree_construction/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 101b_RDRP_phylum_tree.smk
# Run after 03_RDRP_identification.smk (Steps 20-22), 03b, and 03c are
# complete. Three tree-building modes (full_length/palm_core/palm_extended),
# each independently toggleable via config rdrp_tree.{full_length,palm_core,
# palm_extended} -- a single `rule all` builds whichever mode(s) are
# currently enabled, so a bare `snakemake` run with no target builds all of
# them together in one pass.
# ============================================================

# dry-run
snakemake --snakefile workflow/101b_RDRP_phylum_tree.smk --use-conda --cores 1 -n -p

# unlock
snakemake --snakefile workflow/101b_RDRP_phylum_tree.smk --unlock

# formal run
mkdir -p log/101b_RDRP_phylum_tree
nohup snakemake --snakefile workflow/101b_RDRP_phylum_tree.smk \
    --executor slurm \
    --use-conda \
    --jobs 60 \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 10 \
    --keep-going\
    > log/101b_RDRP_phylum_tree/snakemake.log 2>&1 &
echo "PID: $!"



# ============================================================
# 00b_build_diamond_nr_db.smk
# Run on compute node AFTER 00_setup_RNA_database.smk has completed
# (download_blast_nr_db + download_diamond_nr_taxdump need internet;
#  these two rules only read already-downloaded files)
# Step 1: generate_nr_fasta_and_taxid_map  — memory partition, ~12-24 h
# Step 2: build_diamond_nr_db              — compute partition, ~2-4 h
# ============================================================

# dry-run
snakemake --snakefile workflow/00b_build_diamond_nr_db.smk --cores 2 --use-conda -n -p --rerun-triggers input

# build conda envs only
snakemake --snakefile workflow/00b_build_diamond_nr_db.smk --sdm conda --conda-create-envs-only

# unlock
snakemake --snakefile workflow/00b_build_diamond_nr_db.smk --unlock

# formal run (submits two sequential SLURM jobs automatically via snakemake-executor-slurm)
mkdir -p log/00b_build_diamond_nr_db
nohup snakemake --snakefile workflow/00b_build_diamond_nr_db.smk \
    --executor slurm \
    --jobs 2 \
    --use-conda \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 120 \
    > log/00b_build_diamond_nr_db/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# 00_alphafold3.smk
# Singularity/apptainer-based AlphaFold3 install + database download.
# build_alphafold3_sif needs apptainer --fakeroot: run on a compute node via
# an interactive allocation, not the login node (see workflow/00_alphafold3.smk
# for the fallback if fakeroot isn't enabled on this account).
# download_alphafold3_dbs is a long single download (~252 GB / ~630 GB
# uncompressed) -- run it as its own sbatch job, don't do it interactively.
# Model weights (af3.bin.zst) are gated and must be downloaded by hand into
# config["alphafold3"]["model_dir"] -- not part of this workflow.
# ============================================================

# dry-run
snakemake --snakefile workflow/00_alphafold3.smk --cores 1 -n -p

# unlock
snakemake --snakefile workflow/00_alphafold3.smk --unlock


# then, inside that shell:
snakemake --snakefile workflow/00_alphafold3.smk --cores 4  --rerun-triggers input -n -p 

# formal run — db download as its own long sbatch job
mkdir -p log/00_alphafold3
nohup snakemake --snakefile workflow/00_alphafold3.smk \
    --executor slurm \
    --jobs 1 \
    --retries 0 \
    --printshellcmds \
    --slurm-no-account \
    --rerun-triggers input \
    --latency-wait 60 \
    > log/00_alphafold3/snakemake.log 2>&1 &
echo "PID: $!"


# ============================================================
# pipeline figures
# ============================================================
snakemake --snakefile workflow/00_setup_RNA_database.smk           --cores 1 --rulegraph | dot -Tsvg > docs/dag_00_setup_RNA_database.svg
snakemake --snakefile workflow/01_bacteria_genome_assembly.smk     --cores 1 --rulegraph | dot -Tsvg > docs/dag_01_bacteria_genome_assembly.svg
snakemake --snakefile workflow/02_RNA_virus_assembly.smk            --cores 1 --rulegraph | dot -Tsvg > docs/dag_02_RNA_virus_assembly.svg
snakemake --snakefile workflow/02_RNA_virus_assembly_sortmerna.smk  --cores 1 --rulegraph | dot -Tsvg > docs/dag_02_RNA_virus_assembly_sortmerna.svg
snakemake --snakefile workflow/03_RDRP_identification.smk           --cores 1 --rulegraph | dot -Tsvg > docs/dag_03_RDRP_identification.svg
snakemake --snakefile workflow/04_modified_esvirtue.smk             --cores 1 --rulegraph | dot -Tsvg > docs/dag_04_modified_esvirtue.svg

snakemake --snakefile  workflow/03b_esvirtu_rdrp_identification.smk           --cores 1 --rulegraph | dot -Tsvg > docs/03b_esvirtu_rdrp_identification.svg



snakemake --snakefile workflow/05_esvirtue_correlation_analysis.smk --cores 1 --rulegraph | dot -Tsvg > docs/dag_05_esvirtue_correlation_analysis.svg
snakemake --snakefile workflow/101_kraken_check.smk                 --cores 1 --rulegraph | dot -Tsvg > docs/dag_101_kraken_check.svg
