# TET2-CH_CRC-metastasis-scripts

Analysis code for:

**"TET2-mutant clonal hematopoiesis prevents T-cell exhaustion and suppresses cancer metastasis"**

---

## Scripts

### 1. `scRNAseq_complete_pipeline.R`
Single-cell RNA-seq analysis of tumor-infiltrating leukocytes (Vav1-Cre*Tet2*KO vs *Tet2*fl/fl).

| Part | Description | Figure |
|------|-------------|--------|
| 1–2 | Integration (Seurat), dimensionality reduction, whole-leukocyte clustering | — |
| 2B–2C | FindAllMarkers, cluster naming, UMAP/pie/heatmap/violin | Ext. Data Fig. 6f |
| 3 | T cell subset extraction, re-clustering, naming | Ext. Data Fig. 6j |
| 4 | CD8 T cell extraction, re-clustering (T1–T9), heatmap | Ext. Data Fig. 7b |
| 5 | Cell-number equalization (WT = KO) | — |
| 6 | Differential gene expression (volcano plot) | — |
| 7 | AUCell gene signature scoring | — |
| 8 | Exhaustion module score (AddModuleScore) | — |
| 9 | Monocle3 pseudotime trajectory | — |
| 10 | TCR repertoire analysis (scRepertoire) + clonal figures | — |

**R packages required:**
`Seurat`, `SeuratObject`, `dplyr`, `ggplot2`, `patchwork`, `cowplot`, `RColorBrewer`, `ComplexHeatmap`, `circlize`, `openxlsx`, `scRepertoire`, `ggalluvial`, `ggbeeswarm`, `ggpubr`, `viridis`, `tidyr`, `tibble`, `EnhancedVolcano`, `monocle3`, `SingleCellExperiment`, `AUCell`, `babelgene`, `data.table`, `rstatix`

---

### 2. `WGBS_complete_pipeline.R`
Whole-genome bisulfite sequencing analysis — genome-wide DMR calling and annotation (*Tet2*KO vs WT PD-1⁺CD8⁺ TILs).

| Step | Description |
|------|-------------|
| 1–2 | Load Bismark coverage files, QC, DMR calling (methylKit) |
| 3 | Coordinate shift adjustment (NNNNN-trimmed mm39 genome) |
| 4–5 | Build annotation object (exon, intron, promoter, TSS, enhancer, CTCF, UTR, CDS) |
| 6–7 | Annotate DMRs (annotatePeakInBatch + BioMart), gene symbol propagation |
| 8 | Export coordinate-adjusted CSVs |
| 9–10 | Overview figures (stacked bar, pie, volcano, chromosome distribution, violin, lollipop) |
| 11 | Venn diagram — methylation × expression overlap |

**R packages required:**
`methylKit`, `data.table`, `rtracklayer`, `GenomicRanges`, `GenomicFeatures`, `txdbmaker`, `ChIPpeakAnno`, `biomaRt`, `dplyr`, `ggplot2`, `scales`, `patchwork`, `VennDiagram`, `grid`

---

### 3. `Tox_WGBS_complete_pipeline.R`
Focused analysis of the *Tox* locus — CpG-level methylation, DMR characterization, and publication figures.

| Step | Description |
|------|-------------|
| 1–2 | Extract *Tox* CpG sites, filter, normalize, differential methylation (CpG + tile) |
| 3–5 | Annotate *Tox* DMRs, gene symbol propagation, filter biologically relevant DMRs |
| 6–7 | Coordinate adjustment, summary statistics (manuscript numbers) |
| 8 | Fig. — CpG methylation difference track + gene structure overlay |
| 9 | Fig. — DMR scatter plot (spatial landscape) |
| 10 | Fig. — Methylation heatmap (significant CpGs) |
| 11 | Fig. — CpG tracks (WT vs KO individual samples) |
| 12 | Supplementary figures (volcano, DMR characteristics, regulatory element analysis) |

**R packages required:**
`methylKit`, `rtracklayer`, `GenomicRanges`, `ChIPpeakAnno`, `biomaRt`, `dplyr`, `tidyr`, `ggplot2`, `scales`, `patchwork`, `cowplot`, `gridExtra`, `pheatmap`, `viridis`, `ggpubr`, `ggrepel`, `RColorBrewer`, `ComplexHeatmap`, `circlize`

---

## How to use

1. Open each script and edit the **`GLOBAL PATHS`** section at the top — this is the only section that needs to be changed when running on a different machine.
2. Run scripts in this order:
   - `WGBS_complete_pipeline.R` → generates the annotation object (`.rds`) required by `Tox_WGBS_complete_pipeline.R`
   - `Tox_WGBS_complete_pipeline.R`
   - `scRNAseq_complete_pipeline.R` (independent)
3. Intermediate objects are saved as `.rds` checkpoints — subsequent runs will reload from disk automatically.

## Data availability

Raw data are deposited at GEO (accession to be added upon publication).

## Citation

To be added upon publication.
