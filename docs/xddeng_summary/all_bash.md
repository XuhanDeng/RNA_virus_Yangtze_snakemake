Your 48 MT contigs
    │
    ├─ discard < 1000 nt
    ├─ discard rRNA genes (Barrnap or Infernal vs Rfam)
    ├─ dereplicate at 99% (mmseqs easy-linclust)
    │
    ▼
Round 1: MMseqs2 nuc-nuc vs DNAome
    -s 1, --min-seq-id 0.70, --min-aln-len 100, -e 1e-6
    → discard hits, pass survivors forward
    │
    ▼
Round 2: Diamond blastx vs DNAome predicted ORFs
    -k 1, --evalue 1e-5, --query-cover 50
    → discard hits, pass survivors forward
    │
    ▼
Round 3: BLASTn vs DNAome (megablast first, then lower word size)
    -perc_identity 65, -evalue 1e-6
    → discard hits, pass survivors forward
    │
    ▼
Surviving contigs = candidate RNA viruses
    │
    ▼
Secondary filtering (IMG/VR, NR taxonomy, Pfam annotation)