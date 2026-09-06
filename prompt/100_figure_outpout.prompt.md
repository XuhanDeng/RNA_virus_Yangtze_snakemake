# Figure 1
## Network Diagram based on the Spearman Correlation

---

### Input

**Known-unknown pairs (two levels):**

| Level | Directory |
|-------|-----------|
| subspecies | `result/05_esvirtu_correlation_analysis/4_spearman_analysis_RdRp_annotated/subspecies/known_unknown_pair/` |
| species    | `result/05_esvirtu_correlation_analysis/4_spearman_analysis_RdRp_annotated/species/known_unknown_pair/` |

Files: `r0.6.annotated.tsv`, `r0.7.annotated.tsv`, `r0.8.annotated.tsv`, `r0.9.annotated.tsv`

Columns: `seed  target  r  p  r_group  seed_nonzero_samples  target_nonzero_samples  nonzero_intersection  nonzero_union  contig_name  rdrp_category`

- `seed` = known virus (always)
- `target` = unknown contig (always)

Example row:
```
t__Avon-Heathcote Estuary associated circular virus 17  WuRiver_D_3_0000002987  0.904  1.25e-18  0.9  14  12  12  14  WuRiver_D_3_0000002987  not_in_list
```

**Known-known pairs (two levels):**

| Level | Directory |
|-------|-----------|
| subspecies | `result/05_esvirtu_correlation_analysis/3_spearman_analysis/subspecies/known_known_pair/` |
| species    | `result/05_esvirtu_correlation_analysis/3_spearman_analysis/species/known_known_pair/` |

Files: `r0.6.tsv`, `r0.7.tsv`, `r0.8.tsv`, `r0.9.tsv`

Columns: `seed  target  r  p  r_group  seed_nonzero_samples  target_nonzero_samples  nonzero_intersection  nonzero_union`

---

### Node design

| Node type | Color |
|-----------|-------|
| Known virus | one distinct color (e.g. steelblue) |
| Unknown — `High-confidence_RdRp` | distinct color |
| Unknown — `lucaprot` | distinct color |
| Unknown — `RDRPCatch_only` | distinct color |
| Unknown — `RDRPCatch_palmscan` | distinct color |
| Unknown — `not_in_list` | grey |

**Node size** = `nonzero_samples` count for that virus/contig (larger node = detected in more samples).
- For known viruses: use `seed_nonzero_samples` from any row where it appears as seed.
- For unknown contigs: use `target_nonzero_samples`.

---

### Edge design

- Each edge = one correlated pair (from known-unknown or known-known files)
- **Edge width** = proportional to r value
- **Edge color** = r_group threshold bin:
  - 0.6 → color A
  - 0.7 → color B
  - 0.8 → color C
  - 0.9 → color D

---

### Output figures

**16 static figures** (PNG) + **16 interactive figures** (HTML via pyvis):

For each of **2 levels** (subspecies / species) × **4 cumulative thresholds** × **2 label versions**:

| Threshold | Edges included |
|-----------|---------------|
| r ≥ 0.6 | r_group 0.6, 0.7, 0.8, 0.9 |
| r ≥ 0.7 | r_group 0.7, 0.8, 0.9 |
| r ≥ 0.8 | r_group 0.8, 0.9 |
| r ≥ 0.9 | r_group 0.9 only |

Each figure has **two versions**:
- With node labels
- Without node labels

Interactive HTML version (pyvis): hover over a node to see virus name / contig ID / nonzero_samples count.

---

### Output directory structure

```
result/06_network_figure/
  subspecies/
    r0.6_with_label.png
    r0.6_no_label.png
    r0.6_with_label.html
    r0.6_no_label.html
    r0.7_with_label.png
    ...
  species/
    r0.6_with_label.png
    ...
```

---

### Layout

Use **force-directed layout** for all figures:
- Nodes repel each other; edges act as springs pulling correlated nodes together
- Result: viruses/contigs that share many correlations naturally cluster visually without manual arrangement
- For static PNG (`networkx`): use `kamada_kawai_layout` (or `spring_layout` for large graphs)
- For interactive HTML (`pyvis`): force-directed is the default; enable physics so the user can drag and re-settle nodes

---

### Resource config

Use `small_job` config block (lightweight Python script, no GPU needed).
