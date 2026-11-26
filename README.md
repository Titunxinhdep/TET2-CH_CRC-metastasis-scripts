# TET2-CH_CRC-metastasis-scripts

Analysis code for:

**"TET2-mutant clonal hematopoiesis prevents T-cell exhaustion and suppresses cancer metastasis"**

---

## System requirements

**Operating system tested**

- macOS Ventura 13.7 (Intel x86_64)

**R environment**

- R version 4.4.1 (2024-06-14)
- Package versions used for the analyses reported in the manuscript. Where a
  package has since been updated in our environment, the version listed here is
  the one the published results were generated with:

  - `Seurat` 5.5.0
  - `SeuratObject` 5.4.0
  - `dplyr` 1.2.1
  - `ggplot2` 4.0.3
  - `patchwork` 1.3.2
  - `cowplot` 1.2.0
  - `RColorBrewer` 1.1.3
  - `ComplexHeatmap` 2.22.0
  - `circlize` 0.4.18
  - `openxlsx` 4.2.8.1
  - `scRepertoire` 2.2.1
  - `ggalluvial` 0.12.6
  - `ggbeeswarm` 0.7.3
  - `ggpubr` 0.6.3
  - `viridis` 0.6.5
  - `tidyr` 1.3.2
  - `tibble` 3.3.1
  - `EnhancedVolcano` 1.24.0
  - `monocle3` 1.3.7
  - `SingleCellExperiment` 1.28.1
  - `AUCell` 1.28.0
  - `babelgene` 22.9
  - `data.table` 1.18.4
  - `rstatix` 0.7.3
  - `methylKit` 1.30.0
  - `rtracklayer` 1.66.0
  - `GenomicRanges` 1.58.0
  - `GenomicFeatures` 1.58.0
  - `txdbmaker` 1.2.1
  - `ChIPpeakAnno` 3.40.0
  - `biomaRt` 2.62.1
  - `scales` 1.4.0
  - `VennDiagram` 1.8.2
  - `gridExtra` 2.3
  - `pheatmap` 1.0.13
  - `ggrepel` 0.9.8

**Python environment**

- Python 3.9.12 — `numpy` 2.0.2, `pillow` 11.3.0, `scikit-image` 0.24.0, `scipy` 1.13.1 (see `requirements.txt`)

**Scope of this repository**

This repository contains the **downstream analysis** code only. The upstream
processing of the raw sequencing data is not included here and was not run on a
desktop computer:

- **scRNA-seq.** Raw reads were mapped to the mouse genome (GRCm39) and
  demultiplexed with Cell Ranger (10x Genomics, v8.0.1). `scRNAseq_complete_pipeline.R`
  starts from the resulting count matrices.
- **WGBS.** Raw reads were aligned and methylation calls were extracted with the
  Bismark pipeline. `WGBS_complete_pipeline.R` starts from the resulting Bismark
  coverage files.

Both upstream steps were performed on a high-performance computing cluster, as
the intermediate data for a single run reach the tens of terabytes. They cannot
be reproduced on a desktop computer and are therefore outside the scope of these
scripts. The parameters used are described in the Methods section of the
manuscript.

**Non-standard hardware**

For the code in this repository, no non-standard hardware is required and no GPU
is used. All analyses in this repository were performed on an Intel Xeon W
workstation (14 cores, 2.5 GHz) with 128 GB RAM. Genome-wide DMR calling in
`WGBS_complete_pipeline.R` and dataset integration in `scRNAseq_complete_pipeline.R`
are the most memory-intensive steps. `measure_tissue_area.py` runs comfortably
in 8 GB RAM.

The upstream read mapping described above does require high-performance computing
resources and substantial storage.

---

## Installation guide

**R packages** (from an R session):

```r
install.packages("BiocManager")
BiocManager::install(c("methylKit", "rtracklayer", "GenomicRanges", "GenomicFeatures",
  "txdbmaker", "ChIPpeakAnno", "biomaRt", "ComplexHeatmap", "AUCell",
  "SingleCellExperiment", "EnhancedVolcano"))
install.packages(c("Seurat", "SeuratObject", "dplyr", "ggplot2", "patchwork", "cowplot",
  "RColorBrewer", "circlize", "openxlsx", "ggalluvial", "ggbeeswarm", "ggpubr",
  "viridis", "tidyr", "tibble", "data.table", "rstatix", "scales", "VennDiagram",
  "gridExtra", "pheatmap", "ggrepel", "babelgene"))
```

`monocle3` and `scRepertoire` are installed from GitHub:

```r
remotes::install_github("cole-trapnell-lab/monocle3")
remotes::install_github("ncborcherding/scRepertoire")
```

**Python packages:**

```bash
pip install -r requirements.txt
```

**Typical install time** on a normal desktop computer with a broadband
connection: approximately 60 minutes for the full R environment, most of which is
spent compiling Bioconductor dependencies from source, and under 2 minutes for
the Python environment.

---

## Demo

A demonstration dataset is provided in `demo/` for `measure_tissue_area.py`.

**Contents**

```
demo/
└── control_demo/
    └── 4.tif          one stitched whole-section H&E image of a mouse liver,
                       exported from a Keyence BZ-X microscope with the
                       µm-per-pixel calibration retained in the embedded
                       XML metadata
```

`demo/expected_output/` contains the files the script produces from this image,
for comparison.

**Instructions to run on the demo data**

```bash
python3 measure_tissue_area.py "demo"
```

**Expected output**, written to `demo/results/`:

| File | Content |
|------|---------|
| `area_per_image.csv` | one row per image: calibration, tissue area, hole-filled area, number of tissue pieces |
| `total_by_folder.csv` | per-folder totals and image count |
| `missing_calibration.csv` | images with no calibration in the metadata (empty for this demo) |
| `overlay/control_demo_4.jpg` | the input image with the detected tissue boundary drawn in yellow, for visual confirmation |

For this image the script reports:

```
control_demo: 1 images, total 151.48 mm2
```

and `total_by_folder.csv` contains:

```
folder,n_images,total_tissue_mm2,total_tissue_filled_mm2
control_demo,1,151.479,154.052
```

**Expected run time** for the demo: approximately 3 seconds on the workstation
described above (14-core Intel Xeon W). Run time scales with image size.

**Note on the R pipelines.** The single-cell and bisulfite pipelines operate on
full sequencing datasets and cannot be meaningfully demonstrated on a reduced
input. They are run on the processed data deposited at GEO; see the Data
availability section of the manuscript.

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

### 4. `measure_tissue_area.py`

Automated measurement of total liver section area on stitched whole-section H&E
images exported from a Keyence BZ-X microscope. This area is the denominator
used to compute liver metastatic tumour burden; the numerator, tumour area, was
measured separately by manual delineation in Keyence BZ-X Analyzer.

| Step | Description |
|------|-------------|
| 1 | Read the µm-per-pixel calibration from the XML metadata embedded in each BZ-X `.tif` |
| 2 | Convert RGB to HSV and take the saturation channel |
| 3 | Gaussian smoothing (sigma = 2 px) |
| 4 | Global Otsu threshold on the smoothed saturation channel |
| 5 | Discard connected components smaller than 0.1 mm² as debris |
| 6 | Area = retained pixels × pixel size²; a hole-filled variant is recorded in parallel |
| 7 | Write per-image and per-folder CSVs, plus an overlay JPEG with the detected boundary in yellow for visual inspection |

**Usage**

```bash
python3 measure_tissue_area.py "/path/to/parent_folder"
```

The script is pointed at a **parent folder** containing one sub-folder per animal.
Two naming rules are enforced and are easy to miss:

- Only sub-folders whose name begins with `control` or `KO` are read
  (case-insensitive; set by `FOLDER_PREFIXES` at the top of the script).
- Within each sub-folder, only files whose name is a **plain number** are
  measured — `1.tif`, `2.tif` … `10.tif`. Files such as `1-t.tif`, renamed files,
  and macOS `._` companion files are ignored (set by `IMAGE_NAME`).

```
parent_folder/
├── control 6-2/
│   ├── 1.tif      measured
│   ├── 1-t.tif    ignored
│   └── 2.tif      measured
└── KO1-1/
    └── 1.tif      measured
```

Results are written to `parent_folder/results/`.

Individual images excluded because of a technical fault are listed in the
`EXCLUDE` set at the top of the script, as `(folder name, file name)` pairs, so
that every exclusion is visible in the code rather than applied silently. One
image is excluded in this study, `("control 2-2", "9.tif")`: a large region of
the tissue was torn during sectioning, so the section area measured on that
slide would not reflect the true area of the liver section. The exclusion was
made on the basis of the physical condition of the section, independently of the
measured value.

**Python packages required:**
Tested with Python 3.9.12 — `numpy` 2.0.2, `pillow` 11.3.0, `scikit-image` 0.24.0,
`scipy` 1.13.1 (see `requirements.txt`).

---

## How to use

1. Open each script and edit the **`GLOBAL PATHS`** section at the top — this is
   the only section that needs to be changed when running on a different machine.
   This applies to the R scripts only; `measure_tissue_area.py` takes the folder
   as a command-line argument.
2. Run scripts in this order:
   - `WGBS_complete_pipeline.R` → generates the annotation object (`.rds`) required by `Tox_WGBS_complete_pipeline.R`
   - `Tox_WGBS_complete_pipeline.R`
   - `scRNAseq_complete_pipeline.R` (independent)
   - `measure_tissue_area.py` (independent; run from the command line)
3. Intermediate objects are saved as `.rds` checkpoints — subsequent runs will
   reload from disk automatically.

---

## Reproduction instructions

Each table in the Scripts section above maps a numbered part of the pipeline to
the figure panel it generates. To reproduce a given panel, run the parts of the
corresponding script up to and including that number; the intermediate `.rds`
checkpoints make it unnecessary to re-run earlier parts on subsequent passes.

---

## Data availability

Raw data are deposited at GEO (accession to be added upon publication).

---

## License

This repository is released under the MIT License (OSI-approved). See the
`LICENSE` file for the full text.

---

## Citation

To be added upon publication.
