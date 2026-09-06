# parse_motif_hits.R
# Parse and filter RdRP motif search results from RdRP_motif_search.smk
# Output: motif_hits_filtered.tsv, motif_sequences.faa, motif_order.tsv
#
# Requires:
#   - Utils/basicf.r  (GenericHitsParser, Trim2Core2, CalcPcoverage)
#   - R packages: data.table, Biostrings, stringr
#
# Usage:
#   Rscript parse_motif_hits.R <outdir> <input_faa> <rvmt_dir>
#
# Arguments:
#   outdir    : same outdir as used in config.yaml
#   input_faa : your candidate RdRP protein sequences (.faa)
#   rvmt_dir  : path to RVMT-main (for Utils/basicf.r)

args <- commandArgs(trailingOnly = TRUE)

OUTDIR    <- ifelse(length(args) >= 1, args[1], "output")
INPUT_FAA <- ifelse(length(args) >= 2, args[2], "your_RdRPs.faa")
RVMT_DIR  <- ifelse(length(args) >= 3, args[3],
    "/Users/dengxuhan/Library/CloudStorage/OneDrive-DelftUniversityofTechnology/Yangtze_RNA/Script_related/RVMT-main")

suppressPackageStartupMessages({
    library(data.table)
    library(Biostrings)
    library(stringr)
})

source(file.path(RVMT_DIR, "Utils/basicf.r"))

# ─────────────────────────────────────────────────────────────────────────────
# Column name vectors  (from Utils/basicf.r lines 908–916)
# These are already defined by source() above, listed here for reference:
#
#   hhsearch_cols       13 cols: subject_name pCoverage ali_len pL mismatch
#                                gapOpen q1 q2 p1 p2 Probab evalue score
#   hmmsearh_cols       23 cols: domtblout format
#   psiblast_cols       12 cols: pident q1 q2 p1 p2 qL pL ali_len evalue score subject_name
#   MMseq2_outfmt6_cols 19 cols: subject_name evalue ... q1 q2 qL p1 p2 pL ...
#   DaimondP_cols       14 cols: subject_name evalue ... p1 p2 q1 q2 pL qL ...
#
# GenericHitsParser() looks up these vectors via search_tool argument.
# ─────────────────────────────────────────────────────────────────────────────

# Map Snakemake tool folder names → GenericHitsParser search_tool argument
tool_map <- c(
    hhalign   = "hhsearch",
    hmmsearch = "hmmsearch",
    psiblast  = "psiblast",
    mmseqs    = "mmseqs",
    diamond   = "diamondp"
)

motifs <- 1:4
tools  <- names(tool_map)

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3a  Load and parse all 20 raw TSV files (5 tools × 4 motifs)
# ─────────────────────────────────────────────────────────────────────────────
cat("Parsing search results...\n")

all_hits <- rbindlist(lapply(motifs, function(m) {
    rbindlist(lapply(tools, function(tool) {
        f <- file.path(OUTDIR, "search", tool, paste0("mot.", m, ".tsv"))

        if (!file.exists(f)) {
            warning(sprintf("Missing: %s — skipping", f))
            return(NULL)
        }
        if (file.size(f) == 0) {
            message(sprintf("Empty: %s — skipping", f))
            return(NULL)
        }

        h <- tryCatch(
            GenericHitsParser(
                inpt           = f,
                search_tool    = tool_map[tool],
                input_was      = "RdRp_id",
                Query2Profile  = TRUE,
                breakhdrs      = FALSE,
                reducecols     = TRUE,
                calc_pcoverage = TRUE
            ),
            error = function(e) {
                warning(sprintf("Failed to parse %s: %s", f, e$message))
                NULL
            }
        )
        if (is.null(h) || nrow(h) == 0) return(NULL)

        h$motif_type  <- m
        h$search_tool <- tool
        h
    }), fill = TRUE)
}), fill = TRUE)

cat(sprintf("Total raw hits loaded: %d\n", nrow(all_hits)))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3b  Type coercion
# ─────────────────────────────────────────────────────────────────────────────
all_hits[, evalue   := as.numeric(evalue)]
all_hits[, score    := as.numeric(score)]
all_hits[, q1       := as.integer(q1)]
all_hits[, q2       := as.integer(q2)]
all_hits[, p1       := as.integer(p1)]
all_hits[, p2       := as.integer(p2)]
all_hits[, qL       := as.integer(qL)]
all_hits[, pL       := as.integer(pL)]
all_hits[, ali_len  := abs(as.integer(q2) - as.integer(q1))]

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3c  Filter  (Neri thresholds: Annotations.r lines 907–909)
#   evalue < 1e-5
#   score  > 10
#   ali_len > 10 AA  (alignment length on your sequence)
# ─────────────────────────────────────────────────────────────────────────────
hits_filt <- all_hits[evalue  < 1e-5]
hits_filt <- hits_filt[score  > 10  ]
hits_filt <- hits_filt[ali_len > 10 ]

cat(sprintf("Hits after filtering: %d\n", nrow(hits_filt)))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3d  Cull to best hit per sequence per motif (highest score, across tools)
# ─────────────────────────────────────────────────────────────────────────────
setorder(hits_filt, RdRp_id, motif_type, -score, evalue)
hits_best <- hits_filt[, .SD[1], by = .(RdRp_id, motif_type)]

cat(sprintf("Best hits (one per seq per motif): %d\n", nrow(hits_best)))

out_hits <- file.path(OUTDIR, "motif_hits_filtered.tsv")
fwrite(hits_best, out_hits, sep = "\t")
cat(sprintf("Written: %s\n", out_hits))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3e  Extract motif subsequences with Trim2Core2()
#          (Utils/basicf.r line 875)
#
# Trim2Core2() signature:
#   Trim2Core2(hits_df, faa, input_was, subject_was, Expand_projectionX, Out_faa, out_df)
#
#   hits_df           : filtered hit table with RdRp_id, q1, q2, profile columns
#   faa               : AAStringSet of your full-length sequences
#   input_was         : column in hits_df matching names(faa)    → "RdRp_id"
#   subject_was       : column used to build output seq name     → "profile"
#   Expand_projectionX: 0 = exact q1:q2, no flanking expansion
#   Out_faa = TRUE    : return an AAStringSet
#   out_df  = FALSE   : don't return the data frame (set TRUE if you want both)
#
# Output sequence names: "RdRp_id.profile"  e.g. "seq_001.Vfin.mot.1.22"
# ─────────────────────────────────────────────────────────────────────────────
cat("Extracting motif subsequences...\n")

your_faa <- readAAStringSet(INPUT_FAA)

motif_seqs <- Trim2Core2(
    hits_df           = hits_best,
    faa               = your_faa,
    input_was         = "RdRp_id",
    subject_was       = "profile",
    Expand_projectionX = 0,
    Out_faa           = TRUE,
    out_df            = FALSE
)

out_faa <- file.path(OUTDIR, "motif_sequences.faa")
writeXStringSet(motif_seqs, out_faa)
cat(sprintf("Written: %s  (%d sequences)\n", out_faa, length(motif_seqs)))

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4  Determine motif order per sequence
#         Compare q1 positions across mot.1 (A), mot.2 (B), mot.3 (C), mot.4 (D)
#
# Canonical  : A < B < C < D  (normal RdRP)
# Permuted   : C upstream of A  (C-A-B-D configuration)
# Incomplete : fewer than 4 motifs detected
# ─────────────────────────────────────────────────────────────────────────────
cat("Determining motif order...\n")

# Pivot to wide: one row per sequence, q1 for each motif as separate columns
motif_pos <- dcast(
    hits_best,
    RdRp_id ~ motif_type,
    value.var = "q1",
    fun.aggregate = min  # take earliest start if somehow multiple best hits remain
)
setnames(motif_pos,
    old = c("1", "2", "3", "4"),
    new = c("q1_A", "q1_B", "q1_C", "q1_D")
)

# Count how many of the 4 motifs were detected
motif_pos[, n_motifs_found := rowSums(!is.na(.SD)),
          .SDcols = c("q1_A", "q1_B", "q1_C", "q1_D")]

# Canonical: all four present in A < B < C < D order
motif_pos[, canonical := (
    !is.na(q1_A) & !is.na(q1_B) & !is.na(q1_C) & !is.na(q1_D) &
    q1_A < q1_B  & q1_B < q1_C  & q1_C < q1_D
)]

# Permuted: motif C is upstream of motif A  (C-A-B-D)
motif_pos[, permuted := (
    !is.na(q1_C) & !is.na(q1_A) & q1_C < q1_A
)]

# Incomplete: fewer than 4 motifs found
motif_pos[, incomplete := n_motifs_found < 4]

out_order <- file.path(OUTDIR, "motif_order.tsv")
fwrite(motif_pos, out_order, sep = "\t")
cat(sprintf("Written: %s\n", out_order))

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
cat("\n=== Summary ===\n")
cat(sprintf("  Total sequences with ≥1 motif hit : %d\n", nrow(motif_pos)))
cat(sprintf("  Canonical (A-B-C-D, all 4 motifs) : %d\n", sum(motif_pos$canonical,  na.rm = TRUE)))
cat(sprintf("  Permuted  (C-A-B-D)               : %d\n", sum(motif_pos$permuted,   na.rm = TRUE)))
cat(sprintf("  Incomplete (<4 motifs found)       : %d\n", sum(motif_pos$incomplete, na.rm = TRUE)))
cat("\nOutput files:\n")
cat(sprintf("  %s\n", out_hits))
cat(sprintf("  %s\n", out_faa))
cat(sprintf("  %s\n", out_order))
