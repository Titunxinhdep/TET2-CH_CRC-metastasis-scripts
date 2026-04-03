# ==============================================================================
# scRNA-seq COMPLETE ANALYSIS PIPELINE
# Vav1CreTet2KO vs Tet2fl/fl — CD45+ Tumor-Infiltrating Leukocytes
# ==============================================================================
# PIPELINE OVERVIEW:
#   PART 1:   Integration & whole-leukocyte clustering
#   PART 2A:  PCA → UMAP → clustering (all leukocytes)
#   PART 2B:  FindAllMarkers → name leukocyte clusters
#   PART 2C:  Leukocyte figures (UMAP, pie, heatmap, violin)
#             ↳ Extended Data Fig. 6f — leukocyte feature plots
#   PART 3:   T cell subset — clustering, naming, figures
#             ↳ Extended Data Fig. 6j — T cell canonical feature plots
#   PART 4:   CD8 T cell extraction, re-clustering, naming (T1–T9)
#             ↳ Extended Data Fig. 7b — CD8 heatmap (cluster + condition)
#   PART 5:   Cell-number equalization (WT = KO)
#   PART 6:   CD8 equal — DGE (volcano plot)
#   PART 7:   AUCell gene signature scoring
#   PART 8:   Exhaustion module score (AddModuleScore)
#   PART 9:   Monocle3 trajectory
#   PART 10:  TCR repertoire (scRepertoire) + supplementary clonal figures
# ==============================================================================

# ==============================================================================
# GLOBAL PATHS — edit only this section
# ==============================================================================

BASE        <- "/Users/sakatama/Desktop/tran/scRNAseq_Nhan"
OUT_MAIN    <- file.path(BASE, "Combined_cell")
OUT_CD8     <- file.path(BASE, "Combined_cell/CD8 T cells")
OUT_CD8_EQ  <- file.path(BASE, "Combined_cell/CD8 T cells/CD8_equal")
OUT_TRAJ    <- file.path(BASE, "Combined_cell/CD8 T cells/CD8_equal/trajectory/new")
OUT_AUCELL  <- file.path(BASE, "AUCell")
OUT_SUPP    <- file.path(BASE, "Combined_cell/supplemetary")

# TCR contig files
TCR_WT_75 <- "/Users/sakatama/Dropbox/manuscrip-pesonal/scRNA/scRNAseq/cell_rangers/TCR/WT_75_TCR/outs/filtered_contig_annotations.csv"
TCR_WT_89 <- "/Users/sakatama/Dropbox/manuscrip-pesonal/scRNA/scRNAseq/cell_rangers/TCR/WT_89_TCR/outs/filtered_contig_annotations.csv"
TCR_KO_04 <- "/Users/sakatama/Dropbox/manuscrip-pesonal/scRNA/scRNAseq/cell_rangers/TCR/KO_04_TCR/outs/filtered_contig_annotations.csv"
TCR_KO_05 <- "/Users/sakatama/Dropbox/manuscrip-pesonal/scRNA/scRNAseq/cell_rangers/TCR/KO_05_TCR/outs/filtered_contig_annotations.csv"

for (d in c(OUT_MAIN, OUT_CD8, OUT_CD8_EQ, OUT_TRAJ, OUT_AUCELL,
            OUT_SUPP, file.path(OUT_SUPP, "T_cell_defined_cluster"),
            file.path(BASE, "Dropbox/manuscrip-pesonal/sup.figs.TCR.clones"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(openxlsx)
  library(scRepertoire)
  library(ggalluvial)
  library(ggbeeswarm)
  library(ggpubr)
  library(viridis)
  library(tidyr)
  library(tibble)
  library(EnhancedVolcano)
  library(monocle3)
  library(SingleCellExperiment)
  library(AUCell)
  library(babelgene)
  library(data.table)
  library(rstatix)
})

# ==============================================================================
# SHARED AESTHETICS
# ==============================================================================

# Leukocyte colors (B cells, Mono/MΦ/DC, G-MDSC/Gran, T cells, NK cells)
leuko_colors <- c(
  "B cells"    = "#E75480",
  "macrophage" = "#4682B4",
  "granulocyte"= "#3CB371",
  "T cells"    = "#BDB76B",
  "NK cells"   = "#9370DB"
)

# T cell subset colors (CD8, CD4, NK-like T, Atxn+CD4-CD8-, γδT17)
tcell_colors <- c(
  "CD8"           = "#D2691E",
  "CD4"           = "#1E90FF",
  "NK-like T cells"   = "#228B22",
  "Atxn+CD4-CD8-" = "#EFA700",
  "gγδT17"        = "#8B0000"
)

# CD8 cluster colors T1–T9
cd8_colors <- c(
  "T1" = "#9B59B6", "T2" = "#99C4F9", "T3" = "#D98880",
  "T4" = "#F39C12", "T5" = "#3373A6", "T6" = "#B795C0",
  "T7" = "#45B39D", "T8" = "#F7DC6F", "T9" = "#999A2F"
)

condition_colors <- c("WT" = "#5A8BB3", "KO" = "#D87182")

feature_color_ramp <- c("grey90","#FFCCCC","#FF9999","#FF6666",
                        "#FF3333","#FF0000","#CC0000","#990000")

theme_pub <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.line  = element_line(color = "black"),
      axis.ticks = element_line(color = "black"),
      axis.ticks.length = unit(3, "mm"),
      axis.text  = element_text(color = "black"),
      plot.title = element_blank(),
      axis.title = element_blank(),
      legend.title = element_blank(),
      legend.text  = element_blank()
    )
}

# Helper: strip Rps/Rpl/mt genes
filter_genes <- function(seurat_obj) {
  all_g <- rownames(seurat_obj)
  no_rps <- dplyr::setdiff(all_g, grep("^Rps", all_g, value = TRUE))
  no_rpl <- dplyr::setdiff(no_rps, grep("^Rpl", all_g, value = TRUE))
  dplyr::setdiff(no_rpl, grep("^mt-",  all_g, value = TRUE))
}

# ==============================================================================
# PART 1: INTEGRATION (Seurat CCA)
# ==============================================================================

cat("\n======  PART 1: INTEGRATION  ======\n")
# NOTE: WT_75, WT_89 = Tet2fl/fl samples
#       VAV_04, VAV_05 = Vav1CreTet2KO samples
#       (VAV_04 ≡ KO_04, VAV_05 ≡ KO_05 — same samples, different naming convention)

if (!file.exists(file.path(OUT_MAIN, "Combined_cell_TCR_BCR_sample_cellname.rds"))) {

  # Add metadata and unique cell names per sample
  for (obj_name in c("WT_75", "WT_89", "VAV_04", "VAV_05")) {
    obj <- get(obj_name)
    cond  <- ifelse(grepl("WT", obj_name), "WT", "KO")
    obj   <- AddMetaData(obj, metadata = cond,     col.name = "Condition")
    obj   <- AddMetaData(obj, metadata = obj_name, col.name = "Sample")
    obj   <- RenameCells(obj, new.names = paste0(obj_name, "_", colnames(obj)))
    assign(obj_name, obj)
  }

  # CCA integration
  anchors      <- FindIntegrationAnchors(object.list = list(WT_75, WT_89, VAV_04, VAV_05), dims = 1:30)
  Combined_cell <- IntegrateData(anchorset = anchors, dims = 1:30)
  saveRDS(Combined_cell, file.path(OUT_MAIN, "Combined_cell_TCR_BCR_sample_cellname.rds"))

} else {
  Combined_cell <- readRDS(file.path(OUT_MAIN, "Combined_cell_TCR_BCR_sample_cellname.rds"))
  cat("Loaded integrated object.\n")
}

# ==============================================================================
# PART 2A: PCA → UMAP → CLUSTERING (whole leukocytes)
# ==============================================================================

cat("\n======  PART 2A: CLUSTERING  ======\n")

if (!file.exists(file.path(OUT_MAIN, "Combined_cell_sample_cellname_30_1.rds"))) {

  DefaultAssay(Combined_cell) <- "integrated"
  genes_use <- filter_genes(Combined_cell)

  Combined_cell <- ScaleData(Combined_cell, features = genes_use, verbose = FALSE)
  Combined_cell <- RunPCA(Combined_cell, npcs = 200, features = genes_use)
  # QC plots
  ElbowPlot(Combined_cell, ndims = 40)
  JackStrawPlot(Combined_cell <- JackStraw(Combined_cell, num.replicate = 20) %>%
                  ScoreJackStraw(dims = 1:20), dims = 1:20)

  Combined_cell <- RunUMAP(Combined_cell, dims = 1:30)
  Combined_cell <- FindNeighbors(Combined_cell, dims = 1:30)
  Combined_cell <- FindClusters(Combined_cell, resolution = 1)
  saveRDS(Combined_cell, file.path(OUT_MAIN, "Combined_cell_sample_cellname_30_1.rds"))

} else {
  Combined_cell <- readRDS(file.path(OUT_MAIN, "Combined_cell_sample_cellname_30_1.rds"))
  cat("Loaded clustered object.\n")
}

# Visualize
DimPlot(Combined_cell, reduction = "umap", label = TRUE)
DimPlot(Combined_cell, reduction = "umap", split.by = "Condition")

# ==============================================================================
# PART 2B: FIND ALL MARKERS → ANNOTATE LEUKOCYTES → REMOVE CD45-NEG CLUSTERS
# ==============================================================================

cat("\n======  PART 2B: MARKER FINDING & ANNOTATION  ======\n")

DefaultAssay(Combined_cell) <- "RNA"
genes_use <- filter_genes(Combined_cell)

# FindAllMarkers (positive only, for cluster ID)
if (!file.exists(file.path(OUT_MAIN, "exportdata_allmarkers_Combined_cell_30_1_OnlyPositive_500.csv"))) {
  Combined_cell.markers <- FindAllMarkers(Combined_cell, features = genes_use,
                                          only.pos = TRUE, min.pct = 0.25,
                                          logfc.threshold = 0.25)
  write.csv(Combined_cell.markers,
            file.path(OUT_MAIN, "exportdata_allmarkers_Combined_cell_30_1_OnlyPositive_500.csv"))
  # Export per-cluster to Excel
  wb <- createWorkbook()
  for (cl in unique(Combined_cell.markers$cluster)) {
    addWorksheet(wb, sheetName = paste("Cluster", cl))
    writeData(wb, sheet = paste("Cluster", cl),
              Combined_cell.markers[Combined_cell.markers$cluster == cl, ])
  }
  saveWorkbook(wb, file.path(OUT_MAIN,
               "exportdata_allmarkers_Combined_cell_30_1_OnlyPositive_500.xlsx"),
               overwrite = TRUE)
}

# Remove CD45-negative clusters (30, 28, 31) and rename
Combined_cell_New <- subset(Combined_cell, idents = c(30, 28, 31), invert = TRUE)
Combined_cell_New[["RNA"]] <- JoinLayers(Combined_cell_New[["RNA"]])

new.cluster.ids <- c("B cells", "macrophage", "granulocyte",
                     "T cells", "macrophage", "macrophage",
                     "macrophage","B cells","B cells","T cells",
                     "NK cells","T cells", "NK cells", "granulocyte",
                     "T cells", "macrophage", "macrophage","macrophage",
                     "T cells", "T cells", "T cells",
                     "granulocyte","B cells","macrophage", "T cells",
                     "macrophage", "B cells","macrophage", "T cells")
names(new.cluster.ids) <- levels(Combined_cell_New$seurat_clusters)
Combined_cell_New <- RenameIdents(Combined_cell_New, new.cluster.ids)

DefaultAssay(Combined_cell_New) <- "RNA"
saveRDS(Combined_cell_New, file.path(OUT_SUPP, "ALL_COMBINED_CELLS FILTERED/Combined_cell_New.rds"))

# ==============================================================================
# PART 2C: VISUALIZATION — whole leukocytes
# ==============================================================================

cat("\n======  PART 2C: LEUKOCYTE FIGURES  ======\n")

# UMAP
DimPlot(Combined_cell_New, reduction = "umap", pt.size = 0.3,
        cols = leuko_colors, label = FALSE) +
  theme_void() + theme(legend.position = "right",
                       legend.text = element_blank(),
                       legend.title = element_blank())

# Pie charts (WT and KO)
prop_data <- data.frame(cluster = Idents(Combined_cell_New),
                        Condition = Combined_cell_New$Condition) %>%
  group_by(Condition, cluster) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(proportion = count / sum(count) * 100)

make_pie <- function(df) {
  ggplot(df, aes(x = "", y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.25) +
    coord_polar("y") +
    scale_fill_manual(values = leuko_colors) +
    theme_void() + theme(legend.position = "none")
}
make_pie(filter(prop_data, Condition == "WT"))
make_pie(filter(prop_data, Condition == "KO"))
write.csv(prop_data, file.path(OUT_SUPP, "T_cell_defined_cluster/wt_percentage.csv"))

# Heatmap — conserved markers per cell type
DefaultAssay(Combined_cell_New) <- "RNA"
genes_use_new <- filter_genes(Combined_cell_New)

all_conserved_markers_leuko <- list()
for (ct in levels(Idents(Combined_cell_New))) {
  cm <- FindConservedMarkers(Combined_cell_New, features = genes_use_new,
                             ident.1 = ct, grouping.var = "Condition",
                             verbose = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
  log2fc_cols <- grep("avg_log2FC", colnames(cm), value = TRUE)
  cm <- cm[rowSums(cm[log2fc_cols] > 0) == length(log2fc_cols), ]
  cm$gene <- rownames(cm); cm$cell_type <- ct
  all_conserved_markers_leuko[[ct]] <- cm
}
write.csv(do.call(rbind, all_conserved_markers_leuko),
          file.path(OUT_SUPP, "all_conserved_markers_combined.csv"))

desired_order_leuko <- c("B cells","macrophage","granulocyte","T cells","NK cells")
names(leuko_colors) <- desired_order_leuko

top_markers_leuko <- lapply(desired_order_leuko, function(ct) {
  m <- all_conserved_markers_leuko[[ct]]
  head(rownames(m[order(m$max_pval), ]), 20)
})
names(top_markers_leuko) <- desired_order_leuko

# Build average expression matrix
orig_genes   <- unlist(top_markers_leuko)
unique_genes <- unique(orig_genes)
avg_mat      <- as.matrix(AverageExpression(Combined_cell_New,
                                             features = unique_genes,
                                             group.by = "ident",
                                             layer = "data")$RNA)[, desired_order_leuko]
avg_scaled   <- t(scale(t(avg_mat)))
avg_scaled[is.nan(avg_scaled)] <- 0

gene_groups <- rep(desired_order_leuko,
                   sapply(top_markers_leuko, length))
names(gene_groups) <- unlist(top_markers_leuko)

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))
row_annot <- rowAnnotation(
  CellType = anno_simple(gene_groups[rownames(avg_scaled)],
                         col = leuko_colors, width = unit(1.5,"mm"),
                         border = FALSE),
  show_annotation_name = FALSE, show_legend = FALSE)
col_annot <- HeatmapAnnotation(
  CellType = anno_simple(colnames(avg_scaled), col = leuko_colors,
                         pch = 16, height = unit(3,"mm")),
  show_annotation_name = FALSE, show_legend = FALSE)

Heatmap(avg_scaled, col = col_fun,
        top_annotation = col_annot, left_annotation = row_annot,
        show_row_names = FALSE, show_column_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE,
        show_heatmap_legend = FALSE, border = TRUE,
        row_gap = unit(1,"mm"), width = unit(5,"cm"), height = unit(10,"cm"),
        row_split = factor(gene_groups[rownames(avg_scaled)],
                           levels = desired_order_leuko),
        row_title_rot = 0, row_title_gp = gpar(fontsize = 0))

# Stacked violin — canonical markers
canonical_markers_leuko <- c(
  "Cd79a","Cd79b","Ms4a1","Cd19",          # B cells
  "Csf1r","Cd68","Mpeg1","Spi1",           # macrophage
  "S100a8","S100a9","Csf3r","Mmp9",        # granulocyte
  "Cd3e","Cd3d","Cd3g","Trac",             # T cells
  "Ncr1","Klrb1c","Klrk1","Klrd1"         # NK cells
)
VlnPlot(Combined_cell_New, features = canonical_markers_leuko,
        group.by = "ident", pt.size = 0, stack = TRUE, flip = TRUE) +
  theme(axis.text = element_blank(), legend.position = "right",
        legend.text = element_blank())

# ==============================================================================


# B. LEUKOCYTE FEATURE PLOTS — Extended Data Fig. 6f
# Genes (3 per cell type, 5 cell types = 15 genes, 3-column grid):
#   B cells:    Cd19, Cd79a, Cd79b
#   Macrophage: Cd68, Apoc2, Mafb
#   G-MDSC/Gran:Cd33, Csf3r, Clec4d
#   T cells:    Cd3e, Trac, Thy1
#   NK cells:   Ncr1, Klrb1c, Klra7
# Object: Combined_cell_New
# ==============================================================================

genes_6f <- c(
  "Cd19",   "Cd79a",   "Cd79b",    # B cells
  "Cd68",   "Apoc2",   "Mafb",     # Macrophage
  "Cd33",   "Csf3r",   "Clec4d",   # G-MDSC/Granulocyte
  "Cd3e",   "Trac",    "Thy1",     # T cells
  "Ncr1",   "Klrb1c",  "Klra7"     # NK cells
)
genes_6f <- intersect(genes_6f, rownames(Combined_cell_New))
DefaultAssay(Combined_cell_New) <- "RNA"

# Shared feature plot theme
fp_theme <- function() {
  theme_void() +
    theme(
      plot.title        = element_text(size = 12, face = "bold.italic",
                                       hjust = 0.5, color = "black",
                                       margin = margin(t = 5, b = 10)),
      legend.position   = "right",
      legend.key.width  = unit(0.3, "cm"),
      legend.key.height = unit(1,   "cm"),
      legend.title      = element_text(size = 9, angle = 90,
                                       hjust = 0.5, vjust = 0.5),
      legend.text       = element_text(size = 8),
      legend.margin     = margin(l = 10),
      plot.margin       = margin(10, 10, 10, 10),
      panel.background  = element_rect(fill = "white", color = NA),
      plot.background   = element_rect(fill = "white", color = NA)
    )
}

fp_color <- c("grey90","#FFCCCC","#FF9999","#FF6666",
              "#FF3333","#FF0000","#CC0000","#990000")

fp_guide <- guides(colour = guide_colourbar(
  title          = "Expression",
  title.position = "right",
  title.theme    = element_text(angle = 90, hjust = 0.5),
  label.position = "left"))

# Individual PDFs
dir.create(file.path(OUT_SUPP, "feature_plots_leukocyte"),
           showWarnings = FALSE, recursive = TRUE)

for (gene in genes_6f) {
  p <- FeaturePlot(Combined_cell_New, features = gene,
                   reduction = "umap", pt.size = 0.3,
                   order = TRUE, cols = fp_color, raster = FALSE) +
    ggtitle(gene) + fp_theme() + fp_guide
  ggsave(
    file.path(OUT_SUPP, "feature_plots_leukocyte",
              paste0(gene, "_featureplot.pdf")),
    p, width = 4, height = 3.5, dpi = 300, device = "pdf")
  cat("Saved:", gene, "\n")
}

# Combined grid — 3 columns × 5 rows = 15 genes
plots_6f <- lapply(genes_6f, function(gene) {
  FeaturePlot(Combined_cell_New, features = gene,
              reduction = "umap", pt.size = 0.2,
              order = TRUE, cols = fp_color, raster = FALSE) +
    ggtitle(gene) + fp_theme() + fp_guide +
    theme(plot.title = element_text(size = 11, face = "bold.italic"))
})

ggsave(
  file.path(OUT_SUPP, "ExtDataFig6f_Leukocyte_FeaturePlots_grid.pdf"),
  plot   = wrap_plots(plots_6f, ncol = 3),
  width  = 12,
  height = 3.5 * ceiling(length(genes_6f) / 3),
  dpi = 300, device = "pdf")
message("✓ Saved: ExtDataFig6f_Leukocyte_FeaturePlots_grid.pdf")

# ==============================================================================


# PART 3: T CELL SUBSET — clustering, naming, figures
# ==============================================================================

cat("\n======  PART 3: T CELL SUBSET  ======\n")

T_cells <- subset(Combined_cell_New, idents = "T cells")
T_cells@active.ident <- T_cells$seurat_clusters
DefaultAssay(T_cells) <- "RNA"

# Conserved markers per seurat cluster (used to name T cell subtypes)
genes_use_t <- filter_genes(T_cells)
all_conserved_markers_tcells <- list()
for (ct in levels(Idents(T_cells))) {
  cm <- tryCatch(
    FindConservedMarkers(T_cells, features = genes_use_t, ident.1 = ct,
                         grouping.var = "Condition", verbose = FALSE,
                         min.pct = 0.25, logfc.threshold = 0.25),
    error = function(e) NULL)
  if (is.null(cm)) next
  log2fc_cols <- grep("avg_log2FC", colnames(cm), value = TRUE)
  cm <- cm[rowSums(cm[log2fc_cols] > 0) == length(log2fc_cols), ]
  cm$gene <- rownames(cm); cm$cell_type <- ct
  all_conserved_markers_tcells[[ct]] <- cm
}
write.csv(do.call(rbind, all_conserved_markers_tcells),
          file.path(OUT_SUPP, "T_cells_conserved_markers_combined.csv"))

# Rename T cell clusters
new.cluster.ids.tcell <- c("CD8","CD8","CD4",
                           "NK-like T cells","CD4","CD4","Atxn+CD4-CD8-",
                           "gγδT17","CD8")
names(new.cluster.ids.tcell) <- levels(T_cells$seurat_clusters)
T_cells <- RenameIdents(T_cells, new.cluster.ids.tcell)

# UMAP split by Condition
DimPlot(T_cells, reduction = "umap", pt.size = 0.3, cols = tcell_colors,
        split.by = "Condition", label = FALSE) +
  theme_void() + theme(legend.position = "right",
                       legend.text = element_blank(), legend.title = element_blank())

# Stacked bar (% per Condition)
df_tcell <- data.frame(cell_type  = T_cells@active.ident,
                       Condition  = T_cells@meta.data$Condition)
pct_tcell <- df_tcell %>%
  count(cell_type, Condition) %>%
  group_by(Condition) %>%
  mutate(percentage = n / sum(n) * 100) %>% ungroup()

ggplot(pct_tcell, aes(x = Condition, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.2, width = 0.6) +
  scale_fill_manual(values = tcell_colors) +
  scale_y_continuous(limits = c(0,100), expand = c(0,0)) +
  theme_pub() + theme(axis.text.x = element_text(color = "white"),
                      axis.text.y = element_text(color = "white"),
                      axis.ticks = element_line(color = "black"))
write.csv(pct_tcell, file.path(OUT_SUPP, "T_cells_cluster_percentage.csv"))

# T cell heatmap (same approach as leukocyte heatmap — skipping repetitive code)
desired_order_tcell <- c("CD8","CD4","NK-like T cells","Atxn+CD4-CD8-","gγδT17")
names(tcell_colors) <- desired_order_tcell

top_markers_tcell <- lapply(desired_order_tcell, function(ct) {
  m <- all_conserved_markers_tcells[[ct]]
  if (is.null(m)) return(character(0))
  head(rownames(m[order(m$max_pval), ]), 20)
})
names(top_markers_tcell) <- desired_order_tcell

# Feature plots — canonical T cell markers
T_canonical_featureplot <- c(
  "Cd8b1","Cd8a","Nkg7","Tnfrsf4","Cd4","Ass1",
  "Fcer1g","Klrb1c","Klra7","Atxn1","Rabgap1l","Tox",
  "Blk","Il17a","Scart1"
)
FeaturePlot(T_cells, features = T_canonical_featureplot,
            reduction = "umap", pt.size = 0.3, ncol = 5,
            cols = feature_color_ramp, order = TRUE) &
  theme_void() & theme(plot.title = element_text(face = "bold.italic"))

# ==============================================================================


# C. T CELL CANONICAL FEATURE PLOTS — Extended Data Fig. 6j
# Genes (3 per subtype, 5 subtypes = 15 genes, 3-column grid):
#   CD8 T:        Cd8b1, Cd8a, Nkg7
#   CD4 T:        Tnfrsf4, Cd4, Ass1
#   NK-like T:    Fcer1g, Klrb1c, Klra7
#   Atxn+CD4-CD8-:Atxn1, Rabgap1l, Tox
#   γδT17:        Il17a, Blk, Scart1
# Object: T_cells (after T cell subtype naming)
# ==============================================================================

genes_6j <- c(
  "Cd8b1",   "Cd8a",     "Nkg7",      # CD8 T
  "Tnfrsf4", "Cd4",      "Ass1",      # CD4 T
  "Fcer1g",  "Klrb1c",   "Klra7",     # NK-like T
  "Atxn1",   "Rabgap1l", "Tox",       # Atxn+CD4-CD8-
  "Il17a",   "Blk",      "Scart1"     # γδT17
)
genes_6j <- intersect(genes_6j, rownames(T_cells))
DefaultAssay(T_cells) <- "RNA"

# Individual PDFs
dir.create(file.path(OUT_SUPP, "feature_plots_Tcell"),
           showWarnings = FALSE, recursive = TRUE)

for (gene in genes_6j) {
  p <- FeaturePlot(T_cells, features = gene,
                   reduction = "umap", pt.size = 0.3,
                   order = TRUE, cols = fp_color, raster = FALSE) +
    ggtitle(gene) + fp_theme() + fp_guide
  ggsave(
    file.path(OUT_SUPP, "feature_plots_Tcell",
              paste0(gene, "_featureplot.pdf")),
    p, width = 4, height = 3.5, dpi = 300, device = "pdf")
  cat("Saved:", gene, "\n")
}

# Combined grid — 3 columns × 5 rows = 15 genes
plots_6j <- lapply(genes_6j, function(gene) {
  FeaturePlot(T_cells, features = gene,
              reduction = "umap", pt.size = 0.2,
              order = TRUE, cols = fp_color, raster = FALSE) +
    ggtitle(gene) + fp_theme() + fp_guide +
    theme(plot.title = element_text(size = 11, face = "bold.italic"))
})

ggsave(
  file.path(OUT_SUPP, "ExtDataFig6j_Tcell_FeaturePlots_grid.pdf"),
  plot   = wrap_plots(plots_6j, ncol = 3),
  width  = 12,
  height = 3.5 * ceiling(length(genes_6j) / 3),
  dpi = 300, device = "pdf")
message("✓ Saved: ExtDataFig6j_Tcell_FeaturePlots_grid.pdf")


# PART 4: CD8 T CELL EXTRACTION & RE-CLUSTERING
# ==============================================================================

cat("\n======  PART 4: CD8 T CELL EXTRACTION  ======\n")

if (!file.exists(file.path(OUT_CD8, "CD8_T_cells_all_sample_cellname_rawcounts.rds"))) {

  # Extract CD8 clusters (seurat_clusters 3, 9, 29 from Combined_cell_New)
  cluster_cells <- subset(Combined_cell_New, idents = c(3, 9, 29))
  cluster_cells[["RNA"]] <- JoinLayers(cluster_cells[["RNA"]])
  raw_counts  <- GetAssayData(cluster_cells, assay = "RNA")
  CD8_T_cells <- CreateSeuratObject(counts = raw_counts)
  CD8_T_cells@meta.data <- cluster_cells@meta.data[Cells(cluster_cells), ]
  saveRDS(CD8_T_cells,
          file.path(OUT_CD8, "CD8_T_cells_all_sample_cellname_rawcounts.rds"))

} else {
  CD8_T_cells <- readRDS(
    file.path(OUT_CD8, "CD8_T_cells_all_sample_cellname_rawcounts.rds"))
  cat("Loaded CD8 raw object.\n")
}

# Re-clustering CD8 T cells
if (!file.exists(file.path(OUT_CD8, "CD8_T_cells_20_1.0.rds"))) {

  DefaultAssay(CD8_T_cells) <- "RNA"
  genes_cd8 <- filter_genes(CD8_T_cells)

  CD8_T_cells <- NormalizeData(CD8_T_cells) %>%
    FindVariableFeatures(nfeatures = 20000) %>%
    ScaleData(features = genes_cd8) %>%
    RunPCA(features = genes_cd8)

  ElbowPlot(CD8_T_cells, ndims = 30)

  CD8_T_cells <- FindNeighbors(CD8_T_cells, dims = 1:20) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(dims = 1:20)

  saveRDS(CD8_T_cells, file.path(OUT_CD8, "CD8_T_cells_20_1.0.rds"))

} else {
  CD8_T_cells <- readRDS(file.path(OUT_CD8, "CD8_T_cells_20_1.0.rds"))
}

DimPlot(CD8_T_cells, reduction = "umap", label = TRUE, pt.size = 1)
DimPlot(CD8_T_cells, reduction = "umap", label = TRUE, split.by = "Condition")

# FindAllMarkers for CD8 T cells
if (!file.exists(file.path(OUT_CD8, "CD8_T_cells_markers_pos.csv"))) {
  genes_cd8 <- filter_genes(CD8_T_cells)
  markers_cd8 <- FindAllMarkers(CD8_T_cells, features = genes_cd8,
                                only.pos = TRUE, min.pct = 0.1,
                                logfc.threshold = 0.25)
  write.csv(markers_cd8, file.path(OUT_CD8, "CD8_T_cells_markers_pos.csv"))
  wb <- createWorkbook()
  for (cl in unique(markers_cd8$cluster)) {
    addWorksheet(wb, sheetName = paste("Cluster", cl))
    writeData(wb, sheet = paste("Cluster", cl),
              markers_cd8[markers_cd8$cluster == cl, ])
  }
  saveWorkbook(wb, file.path(OUT_CD8, "CD8_T_cells_markers_pos.xlsx"),
               overwrite = TRUE)
}

# Name CD8 clusters T1–T9
# (Names assigned after reviewing FindAllMarkers output and dot plots)
# T1 = TSLC, T2 = Ly6c2hiCD8B TSLC, T3 = TPROG,
# T4 = TEFF-1, T5 = TEFF-2, T6 = TINT, T7 = TTRM, T8 = TTEX, T9 = P-TEX
cd8_cluster_names <- c(
  "0"="T1","1"="T2","2"="T3","3"="T4","4"="T5",
  "5"="T6","6"="T7","7"="T8","8"="T9"
)
CD8_T_cells <- RenameIdents(CD8_T_cells, cd8_cluster_names)
CD8_T_cells$Custom_Cluster  <- Idents(CD8_T_cells)  # store in both columns
CD8_T_cells$named_cluster   <- Idents(CD8_T_cells)  # used by scRepertoire

DimPlot(CD8_T_cells, reduction = "umap", cols = cd8_colors,
        label = FALSE, pt.size = 1) +
  theme_void() + theme(legend.position = "none")

# Dot plot — key markers T1-T9
features_dotplot <- c(
  "Il7r","Ccr7","Sell","Lef1","Tcf7","Id3",          # stem-like
  "Ly6c2","Ifngr2","Nr4a3","Nr4a1","Bhlhe40","Ctla4", # transition
  "Ifng","Ccl4","Klrc1","Cxcr6","Id2","Gzmb","Cd160",
  "Ccr5","Gzma","Cxcr3","Klrk1","Fasl","Itga1","Ccr2","Tnf","Eomes","Irf7",
  "Tox","Ikzf3","Pdcd1","Lag3","Tnfrsf9","Itgav","Dusp4","Cx3cr1",
  "Havcr2","Ikzf2","Gzmk","Itgae",                    # exhaustion
  "Cdk1","Plk1","Aurkb","Ccna2","Ccnb2","Birc5","Mki67" # proliferation
)
DotPlot(CD8_T_cells, features = features_dotplot, group.by = "named_cluster",
        dot.scale = 4) +
  RotatedAxis() +
  scale_color_gradientn(colors = c("#2166AC","white","#B2182B")) +
  theme(axis.text.y = element_text(size = 9))

# ==============================================================================


# A. CD8 HEATMAP — Extended Data Fig. 7b
# Left panel: average expression per cluster (T1–T9)
# Right panel: average expression per Condition (WT vs KO)
# Row labels: italic gene names on right side
# ==============================================================================

# ── 1. Gene list ─────────────────────────────────────────────
genes <- c("Havcr2","Klrc1","Srgap3","Fabp5","Entpd1","Csf1","Pdcd1",
           "Tnfsf4","Gzma","Mki67","Lgals3","Tox","Tigit","Fut8",
           "Tnfrsf9","Il2rb","Lag3","Nab1","Slc2a8","Ctla4","Itm2a",
           "Gcnt1","Psmb8","Plscr1","Kcnk5","Tnip3","Bub1","Psma1",
           "Tank","Sub1","Acot7","Nedd9","Hmgb2","Sdf2l1","Etfb",
           "Ezh2","Mcm3","Mis18bp1","Rad21","Psmd9","Slc43a3","Pmf1")
genes <- intersect(genes, rownames(CD8_T_cells_equal))

# ── 2. Colour settings ───────────────────────────────────────
cluster_colors <- c(
  "T1" = "#9B59B6", "T2" = "#99C4F9", "T3" = "#D98880",
  "T4" = "#F39C12", "T5" = "#3373A6", "T6" = "#B795C0",
  "T7" = "#45B39D", "T8" = "#F7DC6F", "T9" = "#999A2F"
)
condition_colors_hmap <- c("WT" = "#A3C4E0", "KO" = "#F6A3B0")

heatmap_col <- colorRamp2(
  seq(-1, 1, length.out = 100),
  colorRampPalette(c("#05606B", "#FFFFFF", "#9F2704"))(100)
)

# ── 3. Average expression matrices ───────────────────────────

# Left panel: by named_cluster
avg_cluster <- AverageExpression(
  CD8_T_cells_equal, assays = "RNA",
  features = genes, group.by = "named_cluster"
)$RNA
col_order_cluster <- intersect(paste0("T", 1:9), colnames(avg_cluster))
avg_cluster <- avg_cluster[, col_order_cluster, drop = FALSE]
mat_cluster <- t(scale(t(avg_cluster)))
mat_cluster[is.na(mat_cluster)] <- 0
mat_cluster <- mat_cluster[match(genes, rownames(mat_cluster)), , drop = FALSE]

# Right panel: by Condition (WT first, then KO)
avg_cond <- AverageExpression(
  CD8_T_cells_equal, assays = "RNA",
  features = genes, group.by = "Condition"
)$RNA
col_order_cond <- intersect(c("WT", "KO"), colnames(avg_cond))
avg_cond <- avg_cond[, col_order_cond, drop = FALSE]
mat_cond <- t(scale(t(avg_cond)))
mat_cond[is.na(mat_cond)] <- 0
mat_cond <- mat_cond[match(genes, rownames(mat_cond)), , drop = FALSE]

# ── 4. Top annotations ───────────────────────────────────────
top_ann_cluster <- HeatmapAnnotation(
  Cluster = col_order_cluster,
  col = list(Cluster = cluster_colors[col_order_cluster]),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Cluster = list(title = "", nrow = 9, direction = "vertical"))
)

top_ann_cond <- HeatmapAnnotation(
  Condition = col_order_cond,
  col = list(Condition = condition_colors_hmap[col_order_cond]),
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    Condition = list(title = "", nrow = 2, direction = "vertical",
                     labels = c("WT" = "Tet2fl/fl", "KO" = "Vav1CreTet2KO")))
)

# ── 5. Build heatmaps ────────────────────────────────────────
common_hmap_args <- list(
  col = heatmap_col,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  show_column_names = FALSE, show_row_names = FALSE,
  heatmap_legend_param = list(
    title = "Gene Z-score", color_bar = "continuous",
    legend_direction = "horizontal", legend_width = unit(4, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 7)),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height,
              gp = gpar(col = "black", lwd = 0.4, fill = NA))
  }
)

hmap_left <- do.call(Heatmap, c(
  list(matrix = mat_cluster, name = "Gene Z-score",
       top_annotation = top_ann_cluster),
  common_hmap_args
))

hmap_right <- do.call(Heatmap, c(
  list(matrix = mat_cond, name = "Condition",
       top_annotation = top_ann_cond,
       show_heatmap_legend = FALSE),
  common_hmap_args[names(common_hmap_args) != "heatmap_legend_param"]
))

# ── 6. Row labels (italic, right side) ───────────────────────
row_labels <- rowAnnotation(
  Genes = anno_mark(
    at        = seq_along(genes),
    labels    = genes,
    labels_gp = gpar(fontsize = 9, fontface = "italic"),
    labels_rot = 0,
    padding   = unit(1, "mm"),
    link_gp   = gpar(col = "white"),
    link_width  = unit(0, "mm"),
    link_height = unit(0, "mm")
  ),
  width = unit(0.3, "cm") + max_text_width(genes, gp = gpar(fontsize = 9))
)

# ── 7. Draw & save ───────────────────────────────────────────
pdf(file.path(OUT_CD8_EQ, "ExtDataFig7b_CD8_Heatmap_Cluster_Condition.pdf"),
    width = 7, height = 10)
draw(
  hmap_left + hmap_right + row_labels,
  heatmap_legend_side    = "bottom",
  annotation_legend_side = "right",
  merge_legends          = FALSE
)
dev.off()
message("✓ Saved: ExtDataFig7b_CD8_Heatmap_Cluster_Condition.pdf")

# ==============================================================================


# PART 5: CELL-NUMBER EQUALIZATION (WT = KO)
# ==============================================================================

cat("\n======  PART 5: CELL-NUMBER EQUALIZATION  ======\n")

if (!file.exists(file.path(OUT_CD8_EQ, "CD8_T_cells_equal_WT_VAV_equal_cell_number_new.rds"))) {

  set.seed(242)
  wt_idx <- which(CD8_T_cells$Condition == "WT")
  ko_idx <- which(CD8_T_cells$Condition == "KO")
  n_sub  <- min(length(wt_idx), length(ko_idx))

  wt_sub <- subset(CD8_T_cells, cells = sample(wt_idx, n_sub))
  ko_sub <- subset(CD8_T_cells, cells = sample(ko_idx, n_sub))
  cat("WT cells:", n_sub, "| KO cells:", n_sub, "\n")

  CD8_T_cells_equal <- merge(wt_sub, y = ko_sub,
                             add.cell.ids = c("WT","KO"), project = "DEG")
  saveRDS(CD8_T_cells_equal,
          file.path(OUT_CD8_EQ, "CD8_T_cells_equal_WT_VAV_equal_cell_number_new.rds"))

} else {
  CD8_T_cells_equal <- readRDS(
    file.path(OUT_CD8_EQ, "CD8_T_cells_equal_WT_VAV_equal_cell_number_new.rds"))
}
cat("Equal CD8 cells —", table(CD8_T_cells_equal$Condition), "\n")

# Re-normalize, scale, PCA, cluster equal object
if (!file.exists(file.path(OUT_CD8_EQ, "CD8_T_cells_equal_15_1.0.rds"))) {

  DefaultAssay(CD8_T_cells_equal) <- "RNA"
  CD8_T_cells_equal[["RNA"]] <- JoinLayers(CD8_T_cells_equal[["RNA"]])
  genes_eq <- filter_genes(CD8_T_cells_equal)

  CD8_T_cells_equal <- NormalizeData(CD8_T_cells_equal) %>%
    FindVariableFeatures(nfeatures = 20000) %>%
    ScaleData(features = genes_eq) %>%
    RunPCA(features = genes_eq)

  ElbowPlot(CD8_T_cells_equal, ndims = 30)

  CD8_T_cells_equal <- FindNeighbors(CD8_T_cells_equal, dims = 1:15) %>%
    FindClusters(resolution = 1) %>%
    RunUMAP(dims = 1:15)

  # Transfer cluster names from CD8_T_cells to equal object
  # (via cell barcode matching; new clustering may differ — name by marker review)
  Idents(CD8_T_cells_equal) <- CD8_T_cells_equal$RNA_snn_res.1

  saveRDS(CD8_T_cells_equal,
          file.path(OUT_CD8_EQ, "CD8_T_cells_equal_15_1.0.rds"))

} else {
  CD8_T_cells_equal <- readRDS(
    file.path(OUT_CD8_EQ, "CD8_T_cells_equal_15_1.0.rds"))
}

# NOTE: after reviewing FindAllMarkers for the equal object,
#       rename clusters to T1–T9 here:
# cd8_eq_names <- c("0"="T1", ...) ; CD8_T_cells_equal <- RenameIdents(...)
# CD8_T_cells_equal$Custom_Cluster <- Idents(CD8_T_cells_equal)
# CD8_T_cells_equal$named_cluster  <- Idents(CD8_T_cells_equal)

# ==============================================================================
# PART 6: DGE ANALYSIS (KO vs WT) — Volcano plot
# ==============================================================================

cat("\n======  PART 6: DGE VOLCANO  ======\n")

CD8_T_cells_equal[["RNA"]] <- JoinLayers(CD8_T_cells_equal[["RNA"]])
Idents(CD8_T_cells_equal) <- CD8_T_cells_equal$Condition

if (!file.exists(file.path(OUT_CD8_EQ, "CD8_T_cells_equal_dge_results.csv"))) {
  dge <- FindMarkers(CD8_T_cells_equal, ident.1 = "KO", ident.2 = "WT",
                     min.pct = 0.1)
  write.csv(dge, file.path(OUT_CD8_EQ, "CD8_T_cells_equal_dge_results.csv"))
} else {
  dge <- read.csv(file.path(OUT_CD8_EQ, "CD8_T_cells_equal_dge_results.csv"),
                  row.names = 1)
}

dge$color <- "Not Significant"
dge$color[dge$p_val_adj < 0.05 & dge$avg_log2FC > 0] <- "Upregulated KO"
dge$color[dge$p_val_adj < 0.05 & dge$avg_log2FC < 0] <- "Upregulated WT"
color_vec <- c("Upregulated KO" = "#D87182",
               "Upregulated WT" = "#5A8BB3",
               "Not Significant" = "grey")[dge$color]

highlight_genes <- c("Il7r","Ccr7","Sell","Lef1","Tcf7","Id3","Id2",
                     "Gzma","Gzmb","Ifng","Tnf","Eomes","Tbx21","Klrg1",
                     "Pdcd1","Lag3","Havcr2","Tigit","Tox","Prdm1",
                     "Tnfrsf9","Gzmk","Itgae","Entpd1")

EnhancedVolcano(dge,
  lab = rownames(dge), x = "avg_log2FC", y = "p_val_adj",
  title = " ", selectLab = highlight_genes,
  pCutoff = 0.05, FCcutoff = FALSE,
  pointSize = 5, labSize = 5, labFace = "bold", labCol = "black",
  boxedLabels = TRUE, colAlpha = 5/6,
  drawConnectors = TRUE, widthConnectors = 0.5, colConnectors = "black",
  colCustom = color_vec, legendPosition = "none") +
  theme_pub() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(face = "bold", size = 14),
        panel.border = element_rect(color = "black", fill = NA))

# ==============================================================================
# PART 7: AUCell — MSigDB gene set signatures
# ==============================================================================

cat("\n======  PART 7: AUCell SCORING  ======\n")

# Helper: load one MSigDB TSV, convert to mouse orthologs
load_msigdb_geneset <- function(tsv_path, col_name) {
  dt <- fread(tsv_path)
  genes_human <- unlist(strsplit(
    dt[STANDARD_NAME == "GENE_SYMBOLS", ..col_name][[1]], ","))
  orthologs(genes = genes_human, species = "mouse") %>% pull(symbol)
}

aucell_rankings <- AUCell_buildRankings(
  CD8_T_cells_equal@assays$RNA$counts, nCores = 1, splitByBlocks = TRUE)

# Function to run one AUCell comparison and add to metadata
run_aucell <- function(seurat_obj, gene_list_named, score_name,
                       dn_key, up_key, sign = 1) {
  # sign = 1 → score = UP - DN; sign = -1 → score = DN - UP
  auc_mtx <- AUCell_calcAUC(gene_list_named, aucell_rankings, nCores = 1)
  seurat_obj@meta.data <- seurat_obj@meta.data %>%
    rownames_to_column("cells") %>%
    inner_join(as.data.frame(t(getAUC(auc_mtx))) %>%
                 rownames_to_column("cells"), by = "cells") %>%
    column_to_rownames("cells") %>%
    mutate(!!score_name := sign * (get(up_key) - get(dn_key)))
  seurat_obj
}

AUCELL_DIR <- file.path(OUT_AUCELL, "GeneSets")

# 1. Exhausted vs Memory
CD8_T_cells_equal <- run_aucell(
  CD8_T_cells_equal,
  list(
    "CD8 T Exhausted vs Memory DN" = load_msigdb_geneset(
      file.path(AUCELL_DIR, "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_DN.v2024.1.Hs.tsv"),
      "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_DN"),
    "CD8 T Exhausted vs Memory UP" = load_msigdb_geneset(
      file.path(AUCELL_DIR, "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP.v2024.1.Hs.tsv"),
      "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP")
  ),
  score_name = "CD8 T Exhausted vs Memory",
  dn_key = "CD8 T Exhausted vs Memory DN",
  up_key = "CD8 T Exhausted vs Memory UP",
  sign   = -1   # DN - UP → positive = more exhausted than memory
)

# 2. Effector vs Exhausted
CD8_T_cells_equal <- run_aucell(
  CD8_T_cells_equal,
  list(
    "EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN" = load_msigdb_geneset(
      file.path(AUCELL_DIR, "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN.v2024.1.Hs.tsv"),
      "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN"),
    "EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP" = load_msigdb_geneset(
      file.path(AUCELL_DIR, "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP.v2024.1.Hs.tsv"),
      "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP")
  ),
  score_name = "Effector_vs_Exhausted",
  dn_key = "EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN",
  up_key = "EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
  sign   = 1
)

# 3. Naive vs Exhausted
CD8_T_cells_equal <- run_aucell(
  CD8_T_cells_equal,
  list(
    "NAIVE_VS_EXHAUSTED_CD8_TCELL_DN" = load_msigdb_geneset(
      file.path(AUCELL_DIR, "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN.v2024.1.Hs.tsv"),
      "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN"),
    "NAIVE_VS_EXHAUSTED_CD8_TCELL_UP" = load_msigdb_geneset(
      file.path(AUCELL_DIR, "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_UP.v2024.1.Hs.tsv"),
      "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_UP")
  ),
  score_name = "Naive_vs_Exhausted",
  dn_key = "NAIVE_VS_EXHAUSTED_CD8_TCELL_DN",
  up_key = "NAIVE_VS_EXHAUSTED_CD8_TCELL_UP",
  sign   = 1
)

# 4. PD1-high vs PD1-low
CD8_T_cells_equal <- run_aucell(
  CD8_T_cells_equal,
  list(
    "PD1HIGH_VS_PD1LOW_CD8_TCELL_DN" = load_msigdb_geneset(
      file.path(AUCELL_DIR, "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN.v2024.1.Hs.tsv"),
      "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN"),
    "PD1HIGH_VS_PD1LOW_CD8_TCELL_UP" = load_msigdb_geneset(
      file.path(AUCELL_DIR, "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP.v2024.1.Hs.tsv"),
      "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP")
  ),
  score_name = "Pd1_high_vs_Pd1_low",
  dn_key = "PD1HIGH_VS_PD1LOW_CD8_TCELL_DN",
  up_key = "PD1HIGH_VS_PD1LOW_CD8_TCELL_UP",
  sign   = 1
)

saveRDS(CD8_T_cells_equal,
        file.path(OUT_AUCELL, "CD8_T_cells_equal_AUCell_all_scores.rds"))

# AUCell boxplot helper
aucell_boxplot <- function(seurat_obj, score_col, ylim_val,
                           ylab_txt = score_col) {
  seurat_obj@meta.data %>%
    ggboxplot(., "Condition", score_col, color = "Condition",
              order = c("WT","KO"), font.tickslab = 15) +
    geom_quasirandom(aes(color = Condition), width = 0.2, size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    stat_compare_means(comparisons = list(c("WT","KO")),
                       size = 6, tip.length = 0.02) +
    coord_cartesian(ylim = ylim_val) +
    scale_color_manual(values = condition_colors) +
    theme_pub() +
    theme(axis.text = element_blank(),
          axis.ticks = element_line(color = "black"),
          panel.border = element_rect(color = "black", fill = NA))
}

aucell_boxplot(CD8_T_cells_equal, "Effector_vs_Exhausted", c(-0.02, 0.15))
aucell_boxplot(CD8_T_cells_equal, "CD8 T Exhausted vs Memory", c(-0.02, 0.20))
aucell_boxplot(CD8_T_cells_equal, "Naive_vs_Exhausted", c(-0.10, 0.20))
aucell_boxplot(CD8_T_cells_equal, "Pd1_high_vs_Pd1_low", c(-0.10, 0.20))

# AUCell FeaturePlot (split by Condition)
FeaturePlot(CD8_T_cells_equal, features = "Effector_vs_Exhausted",
            reduction = "umap", pt.size = 1, split.by = "Condition") &
  theme_void() &
  scale_color_gradientn(colours = rev(brewer.pal(n=11,"RdBu")))

# ==============================================================================
# PART 8: EXHAUSTION MODULE SCORE (AddModuleScore)
# ==============================================================================

cat("\n======  PART 8: EXHAUSTION SCORE  ======\n")

exhaustion_genes_32 <- c(
  "Plscr1","Srgap3","Tox","Tigit","Havcr2","Lag3",
  "Tnfsf4","Tnfrsf9","Tnip3","Il2rb",
  "Ctla4","Pdcd1","Lgals3","Entpd1",
  "Acot7","Kcnk5","Csf1","Klrc1","Psmb8","Gzma","Hmgb2",
  "Tank","Etfb","Sdf2l1","Psma1","Fabp5","Mki67","Slc43a3","Pmf1",
  "Mis18bp1","Bub1","Mcm3","Slc2a8","Gcnt1","Nab1","Sub1",
  "Psmd9","Fut8","Itm2a","Nedd9","Ezh2","Rad21"
)

Exhaustion_score <- AddModuleScore(
  CD8_T_cells_equal, features = list(exhaustion_genes_32),
  assay = "RNA", name = "Exhaustion_score", seed = 123, slot = "data"
)

# FeaturePlot — exhaustion score on UMAP (split by Condition)
split_plots <- FeaturePlot(Exhaustion_score, features = "Exhaustion_score1",
                           reduction = "umap", pt.size = 1,
                           split.by = "Condition", combine = FALSE)
split_plots <- lapply(split_plots, function(p) {
  p + theme_void() +
    scale_color_viridis(option = "magma", direction = -1,
                        name = "Exhaustion Score") +
    labs(title = "Exhaustion Score on UMAP")
})
wrap_plots(split_plots, ncol = 2)

# Boxplot — exhaustion score WT vs KO
pairwise_test <- wilcox_test(Exhaustion_score1 ~ Condition,
                              data = Exhaustion_score@meta.data) %>%
  adjust_pvalue(method = "bonferroni") %>% add_significance()

Exhaustion_score@meta.data %>%
  ggboxplot(., "Condition", "Exhaustion_score1",
            color = "Condition", order = c("WT","KO")) +
  geom_quasirandom(aes(color = Condition), width = 0.2, size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  stat_compare_means(comparisons = list(c("WT","KO")),
                     method = "wilcox.test", size = 6, tip.length = 0.02) +
  coord_cartesian(ylim = c(-0.5, 1)) +
  scale_color_manual(values = condition_colors) +
  theme_pub() +
  theme(axis.text = element_blank(),
        axis.ticks = element_line(color = "black"),
        panel.border = element_rect(color = "black", fill = NA))

# Save Exhaustion_score object (has Exhaustion_score1 in metadata)
CD8_T_cells_equal$Exhaustion_score1 <- Exhaustion_score$Exhaustion_score1
saveRDS(CD8_T_cells_equal, file.path(OUT_CD8_EQ, "CD8_T_cells_equal_with_scores.rds"))

# ==============================================================================
# PART 9: MONOCLE3 TRAJECTORY
# ==============================================================================

cat("\n======  PART 9: MONOCLE3 TRAJECTORY  ======\n")

CD8_T_cells_equal[["RNA"]] <- JoinLayers(CD8_T_cells_equal[["RNA"]])
expr_mat     <- GetAssayData(CD8_T_cells_equal, assay = "RNA", layer = "data")
cell_meta    <- CD8_T_cells_equal@meta.data
cell_meta    <- cell_meta[colnames(expr_mat), ]
cell_meta$cell_type <- CD8_T_cells_equal@active.ident
gene_meta    <- data.frame(gene_short_name = rownames(expr_mat),
                           row.names       = rownames(expr_mat))

cds <- new_cell_data_set(
  expression_data = as(expr_mat, "sparseMatrix"),
  cell_metadata   = cell_meta,
  gene_metadata   = gene_meta
)

if (!file.exists(file.path(OUT_TRAJ, "cds_dim50_res008.rds"))) {

  cds <- preprocess_cds(cds, num_dim = 50)
  plot_pc_variance_explained(cds)

  cds <- reduce_dimension(cds, reduction_method = "UMAP",
                          preprocess_method = "PCA")
  cds <- cluster_cells(cds, reduction_method = "UMAP",
                       k = 20, cluster_method = "leiden",
                       resolution = 0.008, random_seed = 42)

  # Top markers per partition / cell_type / cluster
  for (grp in c("partition","cell_type","cluster")) {
    mkr <- top_markers(cds, group_cells_by = grp,
                       reduction_method = "UMAP",
                       marker_sig_test = TRUE,
                       genes_to_test_per_group = 200,
                       cores = 4, verbose = FALSE)
    mkr_filt <- mkr[!grepl("Rps|Rpl|mt-", mkr$gene_id), ] %>%
      filter(marker_test_q_value < 0.05) %>%
      arrange(cell_group, desc(pseudo_R2))
    write.csv(mkr_filt,
              file.path(OUT_TRAJ, paste0("top_markers_", grp, ".csv")),
              row.names = FALSE)
  }

  # Learn graph and order cells
  cds <- learn_graph(cds, use_partition = FALSE)
  cds <- order_cells(cds)  # interactive — choose root manually

  saveRDS(cds, file.path(OUT_TRAJ, "cds_dim50_res008.rds"))

} else {
  cds <- readRDS(file.path(OUT_TRAJ, "cds_dim50_res008.rds"))
  cat("Loaded Monocle3 CDS.\n")
}

# Plots
plot_cells(cds, color_cells_by = "cell_type",
           label_cell_groups = FALSE, label_leaves = FALSE,
           label_branch_points = FALSE,
           trajectory_graph_segment_size = 1.5,
           trajectory_graph_color = "#FF0000",
           cell_size = 1) +
  scale_color_manual(values = cd8_colors) +
  theme_void() + theme(legend.text = element_text(size = 5))

plot_cells(cds, color_cells_by = "pseudotime",
           label_cell_groups = FALSE, label_leaves = FALSE,
           label_branch_points = FALSE, label_roots = FALSE,
           trajectory_graph_segment_size = 1.5,
           trajectory_graph_color = "lightgreen",
           cell_size = 1) +
  theme_void()

# ==============================================================================
# PART 10: TCR REPERTOIRE (scRepertoire) + SUPPLEMENTARY CLONAL FIGURES
# ==============================================================================

cat("\n======  PART 10: TCR REPERTOIRE  ======\n")

contig_list <- list(
  read.csv(TCR_WT_75),
  read.csv(TCR_WT_89),
  read.csv(TCR_KO_04),
  read.csv(TCR_KO_05)
)

# NOTE: sample order matches sample IDs used in Seurat integration
#       VAV_04 = KO_04, VAV_05 = KO_05
combined_tcr <- combineTCR(
  contig_list,
  samples  = c("WT_75","WT_89","KO_04","KO_05"),
  ID       = c("WT","WT","KO","KO"),
  removeNA = FALSE, removeMulti = FALSE, filterMulti = TRUE
)

# Integrate TCR info into Seurat object
scRep_example <- combineExpression(
  combined_tcr, CD8_T_cells,
  cloneCall  = "genes",
  group.by   = "sample",
  proportion = FALSE,
  cloneSize  = c(Single = 1, Small = 5, Medium = 20,
                 Large = 100, Hyperexpanded = 500)
)

# Ensure Custom_Cluster / named_cluster are set
scRep_example$Custom_Cluster <- scRep_example$named_cluster

# --- Basic repertoire plots ---
clonalHomeostasis(combined_tcr, cloneCall = "genes", group.by = "ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Clone-size UMAP (all + split by Condition)
colorblind_vector <- c("#4B0055","#009B95","#96D84B","#DF9F0E","#DF7B7B")
DimPlot(scRep_example, group.by = "cloneSize",
        cols = rev(colorblind_vector), pt.size = 1.2) +
  theme_void() + theme(legend.key.size = unit(0.5,"cm"))

DimPlot(scRep_example, group.by = "cloneSize", split.by = "Condition",
        cols = rev(colorblind_vector), pt.size = 1.2) +
  theme_void() + theme(legend.position = "none")

# Clonal overlay
clonalOverlay(scRep_example, reduction = "umap", bins = 10,
              facet = "Condition") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank())

# Chord diagram (clone sharing between clusters)
grid_cols <- cd8_colors[paste0("T", 1:9)]
circles <- getCirclize(scRep_example, group.by = "named_cluster",
                       cloneCall = "genes", proportion = TRUE,
                       include.self = FALSE)
names(grid_cols) <- levels(scRep_example$named_cluster)
chordDiagram(circles, self.link = 1, grid.col = grid_cols,
             directional = 1, direction.type = "arrows",
             link.arr.type = "big.arrow",
             annotationTrack = c("grid","name"),
             annotationTrackHeight = mm_h(c(2.1, 1.1)))

# -----------------------------------------------------------------------
# SUPPLEMENTARY CLONAL FIGURES
# (diversity metrics, Lorenz curves, clone-size distribution,
#  hyperexpanded UMAP, alluvial, gene score vs clone size, overlap)
# -----------------------------------------------------------------------

metadata <- scRep_example@meta.data

# Standardize clone-size factor
clone_levels_orig <- c("Hyperexpanded (100 < X <= 500)",
                       "Large (20 < X <= 100)",
                       "Medium (5 < X <= 20)",
                       "Small (1 < X <= 5)",
                       "Single (0 < X <= 1)")
clone_labels_clean <- c("Hyperexpanded (>100)","Large (21-100)",
                        "Medium (6-20)","Small (2-5)","Single (1)")
metadata$Clone_Size_Clean <- factor(metadata$cloneSize,
                                    levels = clone_levels_orig,
                                    labels = clone_labels_clean)
clone_colors <- c("Single (1)"          = "#4B0055",
                  "Small (2-5)"         = "#009B95",
                  "Medium (6-20)"       = "#96D84B",
                  "Large (21-100)"      = "#DF9F0E",
                  "Hyperexpanded (>100)"= "#DF7B7B")
condition_labels <- c("WT" = "Tet2fl/fl", "KO" = "Tet2KO")

# --- Diversity metrics (Shannon, Inverse Simpson, Gini) ---
calculate_diversity_metrics <- function(meta_df) {
  results <- list()
  for (cond in c("WT","KO")) {
    for (cluster in unique(meta_df$named_cluster)) {
      sub <- meta_df %>%
        filter(Condition == cond, named_cluster == cluster, !is.na(CTaa))
      if (nrow(sub) < 5) next
      freq <- table(sub$CTaa)
      p    <- freq / sum(freq)
      n    <- length(p)
      sorted_f <- sort(as.numeric(freq))
      results[[paste(cond, cluster)]] <- data.frame(
        Condition     = cond,
        Cluster       = as.character(cluster),
        Total_Cells   = nrow(sub),
        Shannon_Index = -sum(p * log(p)),
        Inv_Simpson   = 1 / sum(p^2),
        Gini          = (2 * sum((seq_len(n)) * sorted_f)) /
                        (n * sum(sorted_f)) - (n + 1) / n,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, results)
}

diversity_data <- calculate_diversity_metrics(metadata)
diversity_data$Condition <- factor(diversity_data$Condition, levels = c("WT","KO"))

theme_supp <- theme_classic() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, face = "bold", size = 9),
        legend.position = "top", legend.title = element_blank(),
        axis.text.y = element_text(color = "white"),
        axis.ticks = element_line(color = "black"),
        axis.ticks.length = unit(3,"mm"),
        axis.line = element_line(color = "black"),
        axis.title = element_blank(), plot.title = element_blank())

p_shannon <- ggplot(diversity_data, aes(x = Cluster, y = Shannon_Index, fill = Condition)) +
  geom_bar(stat="identity", position="dodge", color="black", linewidth=0.3) +
  scale_fill_manual(values = condition_colors, labels = condition_labels) +
  geom_text(aes(label = Total_Cells), position = position_dodge(0.9),
            vjust = -0.5, size = 2.5) +
  scale_y_continuous(limits = c(0,6), expand = c(0,0)) + theme_supp

p_inv_simpson <- ggplot(diversity_data, aes(x = Cluster, y = Inv_Simpson, fill = Condition)) +
  geom_bar(stat="identity", position="dodge", color="black", linewidth=0.3) +
  scale_fill_manual(values = condition_colors, labels = condition_labels) +
  scale_y_continuous(limits = c(0,250), expand = c(0,0)) + theme_supp

p_gini <- ggplot(diversity_data, aes(x = Cluster, y = Gini, fill = Condition)) +
  geom_bar(stat="identity", position="dodge", color="black", linewidth=0.3) +
  scale_fill_manual(values = condition_colors, labels = condition_labels) +
  scale_y_continuous(limits = c(0,0.6), expand = c(0,0)) + theme_supp

# --- Lorenz curves ---
lorenz_data <- do.call(rbind, lapply(
  setNames(paste0("T", 1:9), paste0("T", 1:9)),
  function(cl) {
    do.call(rbind, lapply(c("WT","KO"), function(cond) {
      sub <- metadata %>%
        filter(Condition == cond, named_cluster == cl, !is.na(CTaa))
      if (nrow(sub) < 5) return(NULL)
      freq  <- sort(as.numeric(table(sub$CTaa)))
      cs    <- cumsum(freq) / sum(freq)
      data.frame(Condition = cond, Cluster = cl,
                 Prop_Clones = seq_along(freq) / length(freq),
                 Cumsum_Cells = cs)
    }))
  }
))

p_lorenz <- ggplot(lorenz_data,
                   aes(x = Prop_Clones, y = Cumsum_Cells, color = Condition)) +
  geom_line(linewidth = 1) +
  geom_abline(linetype = "dashed", color = "grey50") +
  facet_wrap(~Cluster, nrow = 1) +
  scale_color_manual(values = condition_colors, labels = condition_labels) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) +
  theme_supp +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(color = "white"))

# Combined diversity panel
( p_shannon | p_inv_simpson | p_gini ) / p_lorenz +
  plot_annotation(title = "Clonal Diversity Metrics")

write.csv(diversity_data,
          file.path(OUT_CD8_EQ, "Table_Diversity_Metrics.csv"), row.names = FALSE)

# --- Clone-size distribution stacked bar (WT and KO facets) ---
clone_dist <- metadata %>%
  filter(!is.na(Clone_Size_Clean)) %>%
  count(Condition, named_cluster, Clone_Size_Clean) %>%
  group_by(Condition, named_cluster) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup()

ggplot(clone_dist, aes(x = named_cluster, y = pct, fill = Clone_Size_Clean)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  facet_wrap(~Condition, nrow = 2,
             labeller = labeller(Condition = condition_labels)) +
  scale_fill_manual(values = clone_colors) +
  scale_y_continuous(limits = c(0,110), expand = c(0,0)) +
  theme_supp

# --- Alluvial plot (Medium + Large + Hyperexpanded clones) ---
expanded_clones <- metadata %>%
  filter(!is.na(CTaa)) %>%
  count(CTaa, name = "Clone_Size") %>%
  filter(Clone_Size > 5) %>%
  mutate(Size_Category = case_when(
    Clone_Size <= 20  ~ "Medium (6-20)",
    Clone_Size <= 100 ~ "Large (21-100)",
    TRUE              ~ "Hyperexpanded (>100)"
  ))

alluvial_df <- metadata %>%
  filter(CTaa %in% expanded_clones$CTaa) %>%
  left_join(expanded_clones %>% select(CTaa, Size_Category), by = "CTaa") %>%
  group_by(Condition, named_cluster, Size_Category) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(Condition_Label = condition_labels[Condition])

cond_cols_alluvial <- c("WT" = "#4C7D94", "KO" = "#D5737E")

ggplot(alluvial_df,
       aes(axis1 = Condition_Label, axis2 = named_cluster,
           axis3 = Size_Category, y = n)) +
  geom_alluvium(aes(fill = Condition), width = 1/10, alpha = 0.7) +
  geom_stratum(width = 1/10, fill = "white", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = cond_cols_alluvial) +
  theme_minimal() +
  theme(legend.position = "right", legend.title = element_blank(),
        legend.text = element_blank(), axis.text = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank())

write.csv(alluvial_df,
          file.path(getOption("dropbox",
                              "/Users/sakatama/Dropbox/manuscrip-pesonal"),
                    "sup.figs.TCR.clones/alluvial_summary.csv"),
          row.names = FALSE)

# --- Effector vs Exhausted score × Clone size ---
metadata$Effector_vs_Exhausted <- CD8_T_cells_equal@meta.data[
  rownames(metadata), "Effector_vs_Exhausted"]
metadata$Exhaustion_score1     <- CD8_T_cells_equal@meta.data[
  rownames(metadata), "Exhaustion_score1"]

score_df <- metadata %>%
  filter(!is.na(Clone_Size_Clean)) %>%
  pivot_longer(cols = c(Effector_vs_Exhausted, Exhaustion_score1),
               names_to = "Score_Type", values_to = "Score") %>%
  mutate(Score_Type = gsub("1|_", " ", Score_Type))

ggplot(score_df, aes(x = Clone_Size_Clean, y = Score, fill = Condition)) +
  geom_violin(alpha = 0.6, scale = "width",
              position = position_dodge(width = 0.8)) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8,
               position = position_dodge(width = 0.8)) +
  facet_wrap(~Score_Type, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = condition_colors, labels = condition_labels) +
  scale_y_continuous(expand = c(0.05,0)) +
  theme_supp +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))

# --- Clonotype overlap heatmap (Morisita; WT and KO separately) ---
WT_data <- subset(scRep_example, Condition == "WT")
KO_data <- subset(scRep_example, Condition == "KO")

for (nm in c("WT","KO")) {
  d <- if (nm == "WT") WT_data else KO_data
  ov <- as.matrix(clonalOverlap(d, cloneCall = "genes",
                                 group.by = "named_cluster",
                                 method = "morisita",
                                 exportTable = TRUE))
  cnt <- getCirclize(d, group.by = "named_cluster",
                     cloneCall = "genes", proportion = FALSE,
                     include.self = TRUE)
  # Build combined upper/lower triangle matrix
  full_m <- matrix(NA, nrow = nrow(ov), ncol = ncol(ov),
                   dimnames = dimnames(ov))
  full_m[upper.tri(full_m)] <- ov[upper.tri(ov)]
  # Fill lower triangle from count data
  for (i in seq_len(nrow(cnt))) {
    ri <- which(rownames(full_m) == cnt$from[i])
    ci <- which(colnames(full_m) == cnt$to[i])
    if (ri > ci) full_m[ri, ci] <- cnt$value[i]
    else         full_m[ci, ri] <- cnt$value[i]
  }
  long_m <- full_m %>%
    as.data.frame() %>% rownames_to_column("Row") %>%
    pivot_longer(-Row, names_to = "Column", values_to = "Value") %>%
    mutate(Category = ifelse(is.na(Value), "NA",
                             ifelse(as.character(Row) < Column,
                                    "Frequency", "Count")),
           Row = factor(Row, levels = paste0("T", 1:9)))

  p <- ggplot(long_m,
       aes(x = Column, y = Row,
           fill = ifelse(Category=="Frequency" & !is.na(Value) & Value<1,
                         Value, NA))) +
    geom_tile(color = "black") +
    scale_fill_gradientn(colors = c("#A3C8D2","#FFFFFF","#D65F39"),
                         na.value = "grey90") +
    geom_text(aes(label = ifelse(is.na(Value),"NA",round(Value,3))),
              color = "black", size = 4) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.title = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank()) +
    scale_y_discrete(limits = rev(levels(long_m$Row))) +
    labs(title = paste(nm, "— TCR Overlap (upper: Morisita; lower: Counts)"))
  print(p)
}

cat("\n", strrep("=", 70), "\n")
cat("scRNA-seq PIPELINE COMPLETE\n")
cat(strrep("=", 70), "\n")
