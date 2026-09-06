#!/usr/bin/env bash
# Run after generate_nr_fasta_and_taxid_map completes (nr.fasta exists).
# Deletes BLAST db volume files and tarballs (~1 TB) that are no longer needed.
# Safe to run before build_diamond_nr_db — only nr.fasta and acc2tax are needed.

set -euo pipefail

DB_DIR="/scratch/xddeng/yangtze/RNA/my_rna/database/diamond_nr"
NR_FASTA="$DB_DIR/nr.fasta"

if [[ ! -f "$NR_FASTA" ]]; then
    echo "ERROR: $NR_FASTA not found. Run generate_nr_fasta_and_taxid_map first." >&2
    exit 1
fi

echo "nr.fasta found. Removing BLAST db volumes and tarballs..."

cd "$DB_DIR"

rm -f nr.*.tar.gz nr.*.tar.gz.md5
rm -f nr.*.ppd nr.*.phd nr.*.phi nr.*.pin nr.*.pog nr.*.ppi nr.*.pxm nr.*.psi nr.*.ptf nr.*.pto
rm -f nr.pal nr.pdb nr.pos

echo "Done. Remaining files:"
ls -lh "$DB_DIR"
