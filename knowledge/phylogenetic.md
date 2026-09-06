Here is the full conversation section as markdown:

```markdown
## Adding sequences to the Neri tree

Yes exactly — and this is the complete corrected workflow based on what we
confirmed from the actual EPA-ng documentation.

To directly answer your question — yes, you align your sequences to his
alignment, but not by re-running MAFFT or MUSCLE on everything together.
The alignment happens through the HMM profile, which is built from the
77,510-sequence Neri alignment.

The reason this works is the relationship between the three files:

```
neri_alignment.faa  →  defines column positions  →  neri_tree.nwk
       ↓
  hmmbuild learns those column positions
       ↓
  hmmalign puts YOUR sequences into those same columns
       ↓
  EPA-ng: "your sequences are now in the same coordinate system
           as my reference — I can find where they belong on the tree"
```

The column positions are the common language between the alignment and the
tree. The HMM profile is how you get your new sequences speaking that same
language without disturbing the reference.

---

## Inserting 10,000 ICTV sequences — same workflow

Same workflow exactly — hmmbuild → hmmalign → EPA-ng. The only difference
is the size of your query file. But there is an important question first:

**Why are ICTV sequences not already in the Neri tree?**

Neri 2022 included RefSeq viral sequences as part of their input — so most
ICTV classified viruses should already be in the 77,510-tip tree as
reference sequences. Before doing anything, check:

```bash
# Check if your ICTV sequences are already in the Neri tree
grep ">" neri_alignment.faa | sort > neri_names.txt
grep ">" ictv_rdrps.faa     | sort > ictv_names.txt

# Find which ICTV sequences are NOT already in Neri
comm -23 ictv_names.txt neri_names.txt > ictv_not_in_neri.txt
wc -l ictv_not_in_neri.txt
```

Only insert the ones that are genuinely absent. Inserting sequences already
in the reference tree creates duplicates and distorts branch lengths.

### If they are genuinely new — the workflow is identical

```bash
# Step 1 — same HMM profile you already built (no need to rebuild)

# Step 2 — align ICTV sequences into Neri column space
hmmalign --trim \
         --outformat afa \
         neri_rdrp.hmm \
         ictv_rdrps.faa \
         > ictv_aligned.faa

# Step 3 — run EPA-ng on both sets together
cat your_queries_aligned.faa ictv_aligned.faa > all_queries_aligned.faa

epa-ng \
  --tree neri_tree.nwk \
  --ref-msa neri_alignment.faa \
  --query all_queries_aligned.faa \
  --model WAG+G \
  --outdir epa_output/ \
  --redo \
  --threads 16
```

One query file, one EPA-ng run — your novel vOTUs and the ICTV sequences
get placed simultaneously.

### What is the goal?

| Goal | Best approach |
|---|---|
| Place YOUR novel viruses taxonomically | EPA-ng onto Neri tree — ICTV sequences already in Neri are sufficient anchors |
| Show novel viruses alongside ICTV relatives in a figure | Extract the relevant clade, add your sequences + ICTV relatives, build a dedicated IQ-TREE ML tree |
| Build a new comprehensive reference tree | Rebuild from scratch with IQ-TREE using Neri + ICTV + your sequences |

> For most publication purposes, a dedicated IQ-TREE ML tree for the
> specific clade where your novel sequences sit is what reviewers want to
> see. EPA-ng gives you global placement for taxonomy assignment. The paper
> figure for a novel lineage always shows a proper ML tree with your
> sequences as real tips, not probabilistic placements.

---

## Runtime — 10,000 queries is not a bottleneck

The computational concern is in the wrong place. The 77,510 sequences took
a long time in Neri's pipeline for the **alignment step** — not for EPA-ng.

| Step | Neri's cost | Your EPA-ng cost |
|---|---|---|
| Building the MSA | Weeks — required the whole HHsearch/consensus/expansion pipeline | Zero — Neri's alignment is already done |
| hmmalign | Not applicable | Minutes — sequences aligned independently to the profile |
| FastTree / IQ-TREE | Days — 77,510 tips | Zero — Neri's tree is fixed, never rebuilt |
| EPA-ng placement | Not applicable | Fast — queries are independent |

EPA-ng queries are **embarrassingly parallel** — placing sequence #1 has
zero interaction with placing sequence #5,000.

```bash
epa-ng \
  --tree neri_tree.nwk \
  --ref-msa neri_alignment.faa \
  --query all_10k_queries_aligned.faa \
  --model WAG+G \
  --outdir epa_output/ \
  --redo \
  --threads 32

# hmmalign for 10,000 sequences: ~10–20 minutes
# EPA-ng placement:               ~1–4 hours on 32 threads
# Total: well under one HPC job
```

> What actually takes a long time is building a full IQ-TREE ML tree from
> scratch with all 77,510 Neri + 10,000 ICTV + your novel sequences. EPA-ng
> avoids that entirely — the reference tree is frozen.
```