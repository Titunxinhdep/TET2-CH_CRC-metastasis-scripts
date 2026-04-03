# ==============================================================================
# TOX GENE - WGBS COMPLETE ANALYSIS PIPELINE
# TET2KO vs WT CD8+ T cells - Tox locus methylation analysis
# ==============================================================================
# PIPELINE OVERVIEW:
#   STEP 1:  Extract Tox CpG sites from whole-genome methylation object
#   STEP 2:  Filter, normalize, calculate differential methylation (CpG + tile)
#   STEP 3:  Annotate Tox DMRs with complete annotation object (incl. UTR/CDS)
#   STEP 4:  Gene symbol propagation for enhancer/CTCF
#   STEP 5:  Filter biologically relevant DMRs
#   STEP 6:  Coordinate adjustment (shift back to mm39)
#   STEP 7:  Summary statistics (for manuscript claims)
#   STEP 8:  Figure 1 — CpG methylation difference track + gene structure
#   STEP 9:  Figure 2 — DMR scatter plot (spatial landscape)
#   STEP 10: Figure 3 — Methylation heatmap (significant CpGs)
#   STEP 11: Figure 4 — CpG tracks (WT vs KO individual samples)
#   STEP 12: Supplementary figures (volcano, DMR characteristics,
#              regulatory element analysis, WT vs KO scatter)
# ==============================================================================

# ==============================================================================
# GLOBAL PATHS — edit only this section
# ==============================================================================

BASE      <- "/Users/sakatama/Dropbox/Personal/WGBS/Bismark_pipeline"
OUT_TOX   <- file.path(BASE, "methylKit/Tox/manuscript")

# Canonical annotation object (with UTR + CDS) — output of whole-genome pipeline
RDS_ANNOT_OBJ <- file.path(BASE, "methylKit/annotation_object_complete_with_UTR_CDS.rds")

# Whole-genome methylation object (already built in the main pipeline)
RDS_METH_OBJ <- file.path(BASE, "methylKit/meth_object_CD8.rds")

# GTF (coordinate-shifted)
GTF_SHIFTED <- file.path(BASE, "gencode.vM35.shifted.gtf")

# Coordinate shift table
SHIFT_FILE <- file.path(BASE, "coordinate_shifts_report.csv")

# Tox-specific checkpoints
RDS_TOX_RAW    <- file.path(BASE, "methylKit/Tox/Tox_raw_meth_object_CD8.rds")
RDS_TOX_FILT   <- file.path(BASE, "methylKit/Tox/Tox_meth_object_CD8_lo.count_2_hi.perc_99.9.rds")
RDS_TOX_NORM   <- file.path(BASE, "methylKit/Tox/Tox_meth_normalized_before_unite.rds")
RDS_TOX_UNITED <- file.path(BASE, "methylKit/Tox/Tox_meth_unite_CD8.rds")
RDS_TOX_DIFF_T <- file.path(BASE, "methylKit/Tox/Tox_unfiltered_diffmeth_tile.rds")
RDS_TOX_SIG_DMR<- file.path(BASE, "methylKit/Tox/Tox_significant_DMRs_qvalue_0.01_meth_cutoff_20.rds")
RDS_TOX_DIFF_C <- file.path(BASE, "methylKit/Tox/CpG/Tox_unfiltered_diffmeth_CpG.rds")
RDS_TOX_SIG_DMC<- file.path(BASE, "methylKit/Tox/CpG/Tox_significant_DMRs_CpG_qvalue_0.01_meth_cutoff_20.rds")

# Final output CSVs
CSV_TOX_FINAL     <- file.path(OUT_TOX, "Final_annotated_biologically_relevant_DMRs_Tox.csv")
CSV_TOX_FINAL_ADJ <- file.path(OUT_TOX, "Final_annotated_biologically_relevant_DMRs_Tox_coord_adjusted.csv")

dir.create(OUT_TOX, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BASE, "methylKit/Tox/CpG"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(methylKit)
  library(rtracklayer)
  library(GenomicRanges)
  library(ChIPpeakAnno)
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(cowplot)
  library(gridExtra)
  library(pheatmap)
  library(viridis)
  library(ggpubr)
  library(ggrepel)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
})

# ==============================================================================
# SHARED AESTHETICS
# ==============================================================================

colors_condition <- c("WT" = "#2166AC", "Tet2KO" = "#D6604D")

feature_colors <- c(
  "enhancer" = "#E74C3C",
  "ctcf"     = "#3498DB",
  "tss"      = "#27AE60",
  "TSS"      = "#27AE60",
  "exon"     = "#8E44AD",
  "intron"   = "#2C3E50",
  "promoter" = "#FF7F00",
  "CDS"      = "#FAD7A0",
  "3UTR"     = "#F9E79F",
  "5UTR"     = "#A9CCE3"
)

feature_shapes <- c(
  "enhancer" = 22,
  "ctcf"     = 23,
  "tss"      = 24,
  "TSS"      = 24,
  "exon"     = 25,
  "intron"   = 21,
  "promoter" = 21
)

theme_pub <- theme_minimal() +
  theme(
    text               = element_text(size = 11),
    plot.title         = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title         = element_text(size = 11, face = "bold"),
    axis.text          = element_text(size = 10),
    legend.title       = element_text(size = 10, face = "bold"),
    legend.text        = element_text(size = 9),
    panel.grid.major   = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.border       = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.line          = element_line(color = "black", linewidth = 0.5),
    axis.ticks         = element_line(color = "black", linewidth = 0.5),
    strip.text         = element_text(size = 10, face = "bold"),
    strip.background   = element_rect(color = "black", fill = "grey90", linewidth = 0.5)
  )

# ==============================================================================
# STEP 1-2: TOX CpG EXTRACTION AND DIFFERENTIAL METHYLATION
# ==============================================================================

cat("\n======  STEP 1-2: TOX CpG EXTRACTION & DIFF METHYLATION  ======\n")

if (!file.exists(RDS_TOX_SIG_DMR) || !file.exists(RDS_TOX_SIG_DMC)) {

  meth_obj    <- readRDS(RDS_METH_OBJ)
  gtf_shifted <- import(GTF_SHIFTED)

  # Locate Tox gene
  tox_gene_gr <- gtf_shifted[gtf_shifted$gene_name == "Tox" & gtf_shifted$type == "gene"]
  cat("Tox gene:", as.character(seqnames(tox_gene_gr)), start(tox_gene_gr), "-", end(tox_gene_gr), "\n")

  # Select CpGs overlapping Tox locus
  tox_cpg_raw  <- selectByOverlap(meth_obj, tox_gene_gr)
  tox_cpg_filt <- filterByCoverage(tox_cpg_raw, lo.count = 2, hi.perc = 99.9)
  tox_cpg_norm <- normalizeCoverage(tox_cpg_filt)
  tox_cpg_united <- methylKit::unite(tox_cpg_norm, destrand = FALSE)

  cat("Tox CpG sites after unite:", nrow(getData(tox_cpg_united)), "\n")

  # --- Tile-based DMRs (500 bp windows, step 250 bp) ---
  tiles         <- tileMethylCounts(tox_cpg_united, win.size = 500, step.size = 250,
                                    cov.bases = 0, mc.cores = 4)
  diff_tile     <- calculateDiffMeth(tiles)
  sig_dmr_tox   <- getMethylDiff(diff_tile, difference = 20, qvalue = 0.01)
  cat("Significant tile DMRs:", length(sig_dmr_tox), "\n")

  # --- CpG-level DMCs ---
  diff_cpg      <- calculateDiffMeth(tox_cpg_united)
  sig_dmc_tox   <- getMethylDiff(diff_cpg, difference = 20, qvalue = 0.01)
  cat("Significant DMCs:", length(sig_dmc_tox), "\n")

  saveRDS(tox_cpg_raw,    RDS_TOX_RAW)
  saveRDS(tox_cpg_filt,   RDS_TOX_FILT)
  saveRDS(tox_cpg_norm,   RDS_TOX_NORM)
  saveRDS(tox_cpg_united, RDS_TOX_UNITED)
  saveRDS(diff_tile,      RDS_TOX_DIFF_T)
  saveRDS(sig_dmr_tox,    RDS_TOX_SIG_DMR)
  saveRDS(diff_cpg,       RDS_TOX_DIFF_C)
  saveRDS(sig_dmc_tox,    RDS_TOX_SIG_DMC)

  write.table(as.data.frame(sig_dmr_tox),
              file.path(BASE, "methylKit/Tox/Tox_significant_DMRs_qvalue_0.01_meth_cutoff_20.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(as.data.frame(sig_dmc_tox),
              file.path(BASE, "methylKit/Tox/CpG/Tox_significant_DMRs_CpG_qvalue_0.01_meth_cutoff_20.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

} else {
  tox_cpg_united <- readRDS(RDS_TOX_UNITED)
  sig_dmr_tox    <- readRDS(RDS_TOX_SIG_DMR)
  sig_dmc_tox    <- readRDS(RDS_TOX_SIG_DMC)
  cat("Loaded Tox DMR/DMC objects.\n")
  cat("Tile DMRs:", length(sig_dmr_tox), " | CpG DMCs:", length(sig_dmc_tox), "\n")
}

# ==============================================================================
# STEP 3-5: ANNOTATE TOX DMRs + GENE SYMBOL PROPAGATION + FILTER
# ==============================================================================

cat("\n======  STEP 3-5: ANNOTATE TOX DMRs  ======\n")

if (!file.exists(CSV_TOX_FINAL)) {

  complete_annot <- readRDS(RDS_ANNOT_OBJ)
  cat("Annotation features:", paste(names(complete_annot), collapse = ", "), "\n")

  granges_tox <- as(sig_dmr_tox, "GRanges")
  genome(granges_tox) <- "mm39"

  common_seqlevels <- intersect(seqlevels(granges_tox), seqlevels(complete_annot$promoters))
  granges_tox <- keepSeqlevels(granges_tox, common_seqlevels, pruning.mode = "coarse")

  annotate_feature <- function(peaks, feat_gr, label) {
    res <- annotatePeakInBatch(peaks, AnnotationData = feat_gr)
    mcols(res)$feature_type <- label
    res
  }

  combined_annot_tox <- list(
    promoters = annotate_feature(granges_tox, complete_annot$promoters, "promoter"),
    exons     = annotate_feature(granges_tox, complete_annot$exons,     "exon"),
    introns   = annotate_feature(granges_tox, complete_annot$introns,   "intron"),
    TSS       = annotate_feature(granges_tox, complete_annot$TSS,       "TSS"),
    enhancers = annotate_feature(granges_tox, complete_annot$enhancers, "enhancer"),
    ctcf      = annotate_feature(granges_tox, complete_annot$ctcf,      "ctcf"),
    UTR5      = annotate_feature(granges_tox, complete_annot[["5UTR"]], "5UTR"),
    UTR3      = annotate_feature(granges_tox, complete_annot[["3UTR"]], "3UTR"),
    CDS       = annotate_feature(granges_tox, complete_annot$CDS,       "CDS")
  )

  combined_df_tox <- do.call(rbind, lapply(combined_annot_tox, as.data.frame, row.names = NULL))
  cat("Total annotation rows:", nrow(combined_df_tox), "\n")
  cat("Feature distribution:\n"); print(table(combined_df_tox$feature_type))

  # BioMart gene symbol mapping
  cat("Querying BioMart...\n")
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

  combined_df_tox <- combined_df_tox %>%
    mutate(transcript_id_clean = sub("\\.\\d+$", "", feature))

  tx2gene_full <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                   "external_gene_name", "gene_biotype", "description"),
    filters    = "ensembl_transcript_id",
    values     = unique(combined_df_tox$transcript_id_clean),
    mart       = ensembl
  )

  tx2gene_agg <- tx2gene_full %>%
    group_by(ensembl_transcript_id) %>%
    summarise(
      gene_symbols      = paste(unique(external_gene_name), collapse = ";"),
      gene_biotypes     = paste(unique(gene_biotype),       collapse = ";"),
      gene_descriptions = paste(unique(description),        collapse = ";"),
      .groups = "drop"
    ) %>%
    rename(transcript_id_clean = ensembl_transcript_id)

  combined_df_tox_annotated <- combined_df_tox %>%
    left_join(tx2gene_agg, by = "transcript_id_clean") %>%
    mutate(gene_symbols = ifelse(is.na(gene_symbols) | gene_symbols == "",
                                 feature, gene_symbols))

  # Gene symbol propagation to enhancer/CTCF via overlap
  cat("Propagating gene symbols to enhancer/CTCF...\n")
  gr_all <- makeGRangesFromDataFrame(combined_df_tox_annotated,
                                     seqnames.field = "seqnames",
                                     start.field    = "start",
                                     end.field      = "end",
                                     strand.field   = "strand",
                                     keep.extra.columns = TRUE)

  genic_gr      <- gr_all[gr_all$feature_type %in% c("exon","intron","promoter","TSS","CDS","5UTR","3UTR")]
  regulatory_gr <- gr_all[gr_all$feature_type %in% c("enhancer","ctcf")]

  hits <- findOverlaps(regulatory_gr, genic_gr, ignore.strand = TRUE)
  if (length(hits) > 0) {
    mcols(regulatory_gr)[queryHits(hits), c("gene_symbols","gene_biotypes","gene_descriptions")] <-
      mcols(genic_gr)[subjectHits(hits),  c("gene_symbols","gene_biotypes","gene_descriptions")]
    mcols(regulatory_gr)$annotation_source <- "overlap"
  }

  retained_gr <- gr_all[!gr_all$feature_type %in% c("enhancer","ctcf")]
  final_gr_tox <- c(retained_gr, regulatory_gr)
  final_df_tox <- as.data.frame(final_gr_tox) %>% arrange(seqnames, start)

  # Filter biologically relevant
  meaningful_annotations <- c("inside", "overlapStart", "overlapEnd", "includeFeature")
  biologically_relevant_dmrs_tox <- final_df_tox %>%
    filter(insideFeature %in% meaningful_annotations)

  cat("Unique biologically relevant Tox DMRs:",
      length(unique(biologically_relevant_dmrs_tox$peak)), "\n")
  cat("Feature distribution:\n")
  print(table(biologically_relevant_dmrs_tox$feature_type))

  write.csv(biologically_relevant_dmrs_tox, CSV_TOX_FINAL, row.names = FALSE)
  cat("Saved:", CSV_TOX_FINAL, "\n")

} else {
  biologically_relevant_dmrs_tox <- read.csv(CSV_TOX_FINAL)
  cat("Loaded:", CSV_TOX_FINAL, "\n")
  cat("Unique Tox DMRs:", length(unique(biologically_relevant_dmrs_tox$peak)), "\n")
}

# ==============================================================================
# STEP 6: COORDINATE ADJUSTMENT
# ==============================================================================

cat("\n======  STEP 6: COORDINATE ADJUSTMENT  ======\n")

if (!file.exists(CSV_TOX_FINAL_ADJ)) {

  shift_df <- read.csv(SHIFT_FILE)

  tox_adjusted <- biologically_relevant_dmrs_tox %>%
    left_join(shift_df, by = c("seqnames" = "chr")) %>%
    mutate(
      start          = start          + shift,
      end            = end            + shift,
      start_position = start_position + shift,
      end_position   = end_position   + shift
    ) %>%
    dplyr::select(-shift) %>%
    mutate(meth_status = ifelse(meth.diff > 0, "hyper", "hypo"))

  idx <- which(names(tox_adjusted) == "meth.diff")
  tox_adjusted <- tox_adjusted[, c(1:idx, ncol(tox_adjusted),
                                   (idx+1):(ncol(tox_adjusted)-1))]

  write.csv(tox_adjusted, CSV_TOX_FINAL_ADJ, row.names = FALSE)
  cat("Saved:", CSV_TOX_FINAL_ADJ, "\n")

} else {
  tox_adjusted <- read.csv(CSV_TOX_FINAL_ADJ)
  cat("Loaded:", CSV_TOX_FINAL_ADJ, "\n")
}

# ==============================================================================
# STEP 7: SUMMARY STATISTICS (for manuscript)
# ==============================================================================

cat("\n======  STEP 7: SUMMARY STATISTICS  ======\n")

unique_dmrs_tox <- biologically_relevant_dmrs_tox %>%
  distinct(peak, .keep_all = TRUE)

total_unique    <- nrow(unique_dmrs_tox)
hyper_count     <- sum(unique_dmrs_tox$meth.diff > 0)
hypo_count      <- sum(unique_dmrs_tox$meth.diff < 0)
overall_mean    <- round(mean(unique_dmrs_tox$meth.diff), 2)

enh_ctcf_dmrs <- biologically_relevant_dmrs_tox %>%
  filter(feature_type %in% c("enhancer", "ctcf")) %>%
  distinct(peak, .keep_all = TRUE)
enh_ctcf_mean <- round(mean(enh_ctcf_dmrs$meth.diff), 2)
enh_ctcf_max  <- round(max(enh_ctcf_dmrs$meth.diff),  2)

feature_summary_tox <- biologically_relevant_dmrs_tox %>%
  group_by(feature_type) %>%
  summarise(
    unique_dmrs     = length(unique(peak)),
    mean_meth_diff  = round(mean(meth.diff), 2),
    max_meth_diff   = round(max(meth.diff),  2),
    .groups = "drop"
  ) %>%
  mutate(pct_of_total = round(unique_dmrs / total_unique * 100, 1)) %>%
  arrange(desc(unique_dmrs))

cat("Total unique Tox DMRs:", total_unique, "\n")
cat("All hypermethylated:", all(unique_dmrs_tox$meth.diff > 0), "\n")
cat("Overall mean meth.diff:", overall_mean, "%\n")
cat("Enhancer+CTCF mean meth.diff:", enh_ctcf_mean, "%, max:", enh_ctcf_max, "%\n\n")
cat("Feature breakdown:\n"); print(feature_summary_tox)

write.csv(feature_summary_tox, file.path(OUT_TOX, "Tox_feature_summary.csv"), row.names = FALSE)
write.csv(unique_dmrs_tox,     file.path(OUT_TOX, "Tox_unique_DMRs.csv"),     row.names = FALSE)

# ==============================================================================
# HELPER: PREPARE CpG-LEVEL METHYLATION DATA
# ==============================================================================

build_tox_meth_data <- function(tox_cpg_united) {
  tox_data <- getData(tox_cpg_united)
  data.frame(
    chr    = tox_data$chr,
    start  = tox_data$start,
    end    = tox_data$end,
    WT_02  = tox_data$numCs1 / (tox_data$numCs1 + tox_data$numTs1) * 100,
    WT_04  = tox_data$numCs2 / (tox_data$numCs2 + tox_data$numTs2) * 100,
    WT_06  = tox_data$numCs3 / (tox_data$numCs3 + tox_data$numTs3) * 100,
    KO_08  = tox_data$numCs4 / (tox_data$numCs4 + tox_data$numTs4) * 100,
    KO_10  = tox_data$numCs5 / (tox_data$numCs5 + tox_data$numTs5) * 100,
    KO_12  = tox_data$numCs6 / (tox_data$numCs6 + tox_data$numTs6) * 100
  ) %>%
    mutate(
      WT_mean   = (WT_02 + WT_04 + WT_06) / 3,
      KO_mean   = (KO_08 + KO_10 + KO_12) / 3,
      meth_diff = KO_mean - WT_mean,
      position  = start
    )
}

tox_meth_data <- build_tox_meth_data(tox_cpg_united)
x_range <- range(tox_meth_data$position)

# ==============================================================================
# HELPER: BUILD GENE STRUCTURE ANNOTATION FROM ANNOTATION OBJECT
# ==============================================================================

build_tox_gene_annotation <- function(complete_annot, gtf_shifted, x_range) {

  tox_gene <- gtf_shifted[gtf_shifted$gene_name == "Tox" & gtf_shifted$type == "gene"][1]
  tox_chr  <- as.character(seqnames(tox_gene))
  ext_region <- GRanges(seqnames = tox_chr,
                        ranges = IRanges(start = start(tox_gene) - 50000,
                                         end   = end(tox_gene)   + 50000))

  feature_df <- data.frame()
  for (feat_name in names(complete_annot)) {
    feat_gr  <- complete_annot[[feat_name]]
    overlaps <- subsetByOverlaps(feat_gr, ext_region)
    if (length(overlaps) == 0) next
    df <- data.frame(
      chr          = as.character(seqnames(overlaps)),
      start        = start(overlaps),
      end          = end(overlaps),
      feature_type = tolower(feat_name),
      mid_position = (start(overlaps) + end(overlaps)) / 2,
      stringsAsFactors = FALSE
    )
    feature_df <- rbind(feature_df, df)
  }

  # Standardize names
  feature_df <- feature_df %>%
    mutate(feature_type = case_when(
      feature_type == "tss"      ~ "tss",
      feature_type == "exons"    ~ "exon",
      feature_type == "introns"  ~ "intron",
      feature_type == "enhancers"~ "enhancer",
      feature_type == "ctcf"     ~ "ctcf",
      TRUE ~ feature_type
    )) %>%
    filter(feature_type %in% names(feature_colors),
           start <= max(x_range), end >= min(x_range))

  feature_df
}

# Build annotation for plots
if (exists("complete_annot")) {
  gtf_shifted <- if (exists("gtf_shifted")) gtf_shifted else import(GTF_SHIFTED)
  tox_gene_annot <- build_tox_gene_annotation(complete_annot, gtf_shifted, x_range)
} else {
  complete_annot <- readRDS(RDS_ANNOT_OBJ)
  gtf_shifted    <- import(GTF_SHIFTED)
  tox_gene_annot <- build_tox_gene_annotation(complete_annot, gtf_shifted, x_range)
}

# ==============================================================================
# HELPER: CREATE GENE TRACK PANEL (aligned to x_range)
# ==============================================================================

create_gene_track <- function(tox_gene_annot, x_range, height = 0) {

  exons     <- tox_gene_annot %>% filter(feature_type == "exon")   %>% distinct(start, end, .keep_all = TRUE)
  enhancers <- tox_gene_annot %>% filter(feature_type == "enhancer") %>% distinct(start, end, .keep_all = TRUE)
  ctcf_feat <- tox_gene_annot %>% filter(feature_type == "ctcf")   %>% distinct(start, end, .keep_all = TRUE)
  tss_feat  <- tox_gene_annot %>% filter(feature_type %in% c("tss","TSS")) %>% distinct(start, end, .keep_all = TRUE)
  introns   <- tox_gene_annot %>% filter(feature_type == "intron") %>% distinct(start, end, .keep_all = TRUE)

  # Gene backbone line
  p <- ggplot() +
    annotate("segment", x = min(x_range), xend = max(x_range),
             y = height, yend = height, linewidth = 1.5, color = "black")

  # Exons
  if (nrow(exons) > 0)
    p <- p + geom_rect(data = exons,
                       aes(xmin = start, xmax = end,
                           ymin = height - 5, ymax = height + 5),
                       fill = "#8B0000", color = NA, alpha = 0.9, inherit.aes = FALSE)

  # Introns (dashed segments between exons, filtered to non-exon regions)
  if (nrow(introns) > 0 && nrow(exons) > 0) {
    intron_gr <- GRanges(seqnames = introns$chr,
                         ranges = IRanges(introns$start, introns$end))
    exon_gr   <- GRanges(seqnames = exons$chr,
                         ranges = IRanges(exons$start, exons$end))
    non_exon_introns <- as.data.frame(subsetByOverlaps(intron_gr, exon_gr, invert = TRUE))
    if (nrow(non_exon_introns) > 0)
      p <- p + geom_segment(data = non_exon_introns,
                             aes(x = start, xend = end, y = height, yend = height),
                             linewidth = 1.2, color = "gray60",
                             linetype = "dashed", inherit.aes = FALSE)
  }

  # Enhancers (downward triangles above gene line)
  if (nrow(enhancers) > 0)
    p <- p + geom_point(data = enhancers,
                        aes(x = mid_position, y = height + 4),
                        shape = 25, size = 4, fill = "grey50",
                        color = "black", stroke = 0.3, inherit.aes = FALSE)

  # TSS (upward triangles below gene line)
  if (nrow(tss_feat) > 0)
    p <- p + geom_point(data = tss_feat,
                        aes(x = mid_position, y = height - 5.2),
                        shape = 24, size = 4, fill = "grey90",
                        color = "black", stroke = 0.3, inherit.aes = FALSE)

  # CTCF (diamonds on gene line)
  if (nrow(ctcf_feat) > 0)
    p <- p + geom_point(data = ctcf_feat,
                        aes(x = mid_position, y = height + 0.5),
                        shape = 23, size = 4, fill = "grey30",
                        color = "black", stroke = 0.3, inherit.aes = FALSE)

  p <- p +
    coord_cartesian(xlim = x_range, ylim = c(-12, 12)) +
    theme_pub +
    theme(axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          plot.margin  = margin(2, 20, 5, 20)) +
    labs(x = "Genomic position (bp)", y = NULL)

  p
}

gene_track <- create_gene_track(tox_gene_annot, x_range)

# ==============================================================================
# STEP 8: FIGURE 1 — CpG METHYLATION DIFFERENCE TRACK + GENE STRUCTURE
# ==============================================================================

cat("\n======  STEP 8: FIGURE 1 — CpG DIFF TRACK  ======\n")

# Panel A: methylation difference per CpG (line + points)
fig1_panel_a <- tox_meth_data %>%
  ggplot(aes(x = position, y = meth_diff)) +
  geom_line(color = "#e74c3c", linewidth = 0.8) +
  geom_point(color = "#e74c3c", alpha = 0.6, size = 1.2) +
  geom_hline(yintercept = c(-20, 20), linetype = "dashed",
             color = "black", alpha = 0.7, linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "solid",
             color = "black", alpha = 0.5, linewidth = 0.4) +
  coord_cartesian(xlim = x_range, ylim = c(-70, 100)) +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Methylation difference (%)\n(Tet2KO - WT)",
       subtitle = sprintf("All CpG sites at Tox locus (n = %d)", nrow(tox_meth_data)))

# Panel B: gene structure track
fig1_panel_b <- gene_track

# Panel C: boxplot + jitter by feature type
dmr_for_plot <- biologically_relevant_dmrs_tox %>%
  filter(!is.na(feature_type),
         feature_type %in% names(feature_colors)) %>%
  mutate(feature_type = factor(feature_type, levels = names(feature_colors)))

fig1_panel_c <- ggplot(dmr_for_plot, aes(x = feature_type, y = meth.diff)) +
  geom_boxplot(fill = NA, color = "black", outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(width = 0.2, size = 3, alpha = 0.8,
              shape = 21, color = "black", fill = "grey40") +
  geom_hline(yintercept = c(-20, 20), linetype = "dashed",
             color = "black", alpha = 0.7) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Feature type", y = "Methylation difference (%)")

figure_1 <- plot_grid(
  plot_grid(fig1_panel_a, fig1_panel_b, ncol = 1,
            rel_heights = c(0.75, 0.25), align = "v"),
  fig1_panel_c,
  labels   = c("A", "C"),
  ncol     = 2,
  rel_widths = c(4, 1),
  align    = "h"
)

ggsave(file.path(OUT_TOX, "Fig1_Tox_CpG_diff_track.pdf"), figure_1, width = 14, height = 6, dpi = 300)
ggsave(file.path(OUT_TOX, "Fig1_Tox_CpG_diff_track.png"), figure_1, width = 14, height = 6, dpi = 300)
cat("Saved: Fig1_Tox_CpG_diff_track\n")

# ==============================================================================
# STEP 9: FIGURE 2 — SPATIAL LANDSCAPE (DMR scatter by position)
# ==============================================================================

cat("\n======  STEP 9: FIGURE 2 — DMR SPATIAL LANDSCAPE  ======\n")

plot_data_dmr <- biologically_relevant_dmrs_tox %>%
  filter(!is.na(feature_type), feature_type %in% names(feature_colors)) %>%
  mutate(
    position       = (start + end) / 2,
    abs_meth_diff  = abs(meth.diff),
    feature_type   = factor(feature_type, levels = names(feature_colors))
  )

fig2_panel_a <- ggplot(plot_data_dmr,
                        aes(x = position, y = meth.diff)) +
  geom_point(aes(fill = feature_type, shape = feature_type, size = abs_meth_diff),
             color = "black", stroke = 0.5, alpha = 0.9,
             position = position_jitter(width = 3000, height = 1.5)) +
  scale_fill_manual(values  = feature_colors, name = "Feature type") +
  scale_shape_manual(values = feature_shapes, name = "Feature type") +
  scale_size_continuous(name = "|Δmeth| (%)", range = c(2, 7)) +
  geom_hline(yintercept = c(-20, 20), linetype = "dashed",
             color = "black", alpha = 0.6) +
  coord_cartesian(xlim = x_range, ylim = c(-40, 70)) +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  guides(fill  = guide_legend(override.aes = list(size = 5, stroke = 0.5, color = "black")),
         shape = guide_legend(override.aes = list(size = 5, stroke = 0.5, color = "black"))) +
  labs(y = "Methylation difference (%)")

figure_2 <- plot_grid(
  fig2_panel_a, gene_track,
  ncol = 1, rel_heights = c(0.75, 0.25), align = "v",
  labels = c("A", "")
)

ggsave(file.path(OUT_TOX, "Fig2_Tox_DMR_landscape.pdf"), figure_2, width = 12, height = 7, dpi = 300)
ggsave(file.path(OUT_TOX, "Fig2_Tox_DMR_landscape.png"), figure_2, width = 12, height = 7, dpi = 300)
cat("Saved: Fig2_Tox_DMR_landscape\n")

# ==============================================================================
# STEP 10: FIGURE 3 — METHYLATION HEATMAP (significant CpGs)
# ==============================================================================

cat("\n======  STEP 10: FIGURE 3 — METHYLATION HEATMAP  ======\n")

sig_dmc_df    <- as.data.frame(getData(sig_dmc_tox))
sig_positions <- sig_dmc_df$start

meth_matrix_all <- getData(tox_cpg_united) %>%
  as.data.frame() %>%
  mutate(
    WT_02 = numCs1 / (numCs1 + numTs1) * 100,
    WT_04 = numCs2 / (numCs2 + numTs2) * 100,
    WT_06 = numCs3 / (numCs3 + numTs3) * 100,
    KO_08 = numCs4 / (numCs4 + numTs4) * 100,
    KO_10 = numCs5 / (numCs5 + numTs5) * 100,
    KO_12 = numCs6 / (numCs6 + numTs6) * 100
  ) %>%
  filter(start %in% sig_positions, complete.cases(.))

heatmap_matrix           <- as.matrix(meth_matrix_all[, c("WT_02","WT_04","WT_06","KO_08","KO_10","KO_12")])
rownames(heatmap_matrix) <- paste0("CpG_", meth_matrix_all$start)
colnames(heatmap_matrix) <- c("WT-02","WT-04","WT-06","KO-08","KO-10","KO-12")

sample_annot <- data.frame(
  Condition = c(rep("WT",3), rep("Tet2KO",3)),
  row.names = colnames(heatmap_matrix)
)
annot_colors <- list(Condition = c("WT" = "#2166AC", "Tet2KO" = "#D6604D"))

heatmap_colors <- colorRampPalette(
  c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0",
    "#F7F7F7","#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))(100)

pdf(file.path(OUT_TOX, "Fig3_Tox_methylation_heatmap.pdf"), width = 6, height = 8)
pheatmap(
  heatmap_matrix,
  annotation_col          = sample_annot,
  annotation_colors       = annot_colors,
  color                   = heatmap_colors,
  scale                   = "row",
  clustering_distance_rows= "correlation",
  clustering_method       = "ward.D2",
  cluster_rows            = TRUE,
  cluster_cols            = FALSE,
  cutree_rows             = 3,
  show_rownames           = FALSE,
  show_colnames           = FALSE,
  treeheight_row          = 0,
  treeheight_col          = 0,
  main                    = "",
  border_color            = NA,
  legend                  = TRUE,
  annotation_legend       = TRUE,
  annotation_names_col    = FALSE
)
dev.off()
cat("Saved: Fig3_Tox_methylation_heatmap.pdf\n")

# ==============================================================================
# STEP 11: FIGURE 4 — CpG TRACKS (individual samples WT vs KO)
# ==============================================================================

cat("\n======  STEP 11: FIGURE 4 — CpG INDIVIDUAL TRACKS  ======\n")

tox_raw_df <- getData(tox_cpg_united) %>% as.data.frame()
n_samples  <- (ncol(tox_raw_df) - 4) / 3
treatment  <- getTreatment(tox_cpg_united)
sample_ids <- getSampleID(tox_cpg_united)

tracks_list <- lapply(seq_len(n_samples), function(i) {
  data.frame(
    pos          = tox_raw_df$start,
    methylation  = tox_raw_df[[paste0("numCs",i)]] /
                   tox_raw_df[[paste0("coverage",i)]] * 100,
    sample       = sample_ids[i],
    sample_type  = ifelse(treatment[i] == 0, "WT", "KO"),
    stringsAsFactors = FALSE
  )
})
tracks_data <- do.call(rbind, tracks_list)

# Highlight significant DMR windows
sig_dmr_df  <- as.data.frame(getData(sig_dmr_tox))
plot_range  <- range(tracks_data$pos)
sig_in_range <- sig_dmr_df %>%
  filter(start >= plot_range[1], end <= plot_range[2])

fig4_tracks <- ggplot(tracks_data, aes(x = pos, y = methylation, color = sample_type)) +
  geom_point(size = 1.2, alpha = 0.7) +
  facet_wrap(~sample_type, ncol = 1) +
  scale_color_manual(values = c("WT" = "#2166AC", "KO" = "#D6604D")) +
  ylim(0, 100) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Genomic position (bp)", y = "Methylation (%)",
       title = "Tox locus CpG methylation — individual samples")

if (nrow(sig_in_range) > 0) {
  sig_rect <- data.frame(xmin = sig_in_range$start, xmax = sig_in_range$end,
                         ymin = -5, ymax = 105)
  fig4_tracks <- fig4_tracks +
    geom_rect(data = sig_rect,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "grey50", alpha = 0.2, color = "grey50",
              linewidth = 0.5, inherit.aes = FALSE)
}

# LOESS trend lines
group_means <- tracks_data %>%
  group_by(sample_type, pos) %>%
  summarise(mean_meth = mean(methylation, na.rm = TRUE), .groups = "drop")

smooth_lines <- group_means %>%
  group_by(sample_type) %>%
  group_modify(~{
    lo <- tryCatch(loess(mean_meth ~ pos, data = .x, span = 0.3), error = function(e) NULL)
    grid_pos <- seq(min(.x$pos), max(.x$pos), length.out = 800)
    if (is.null(lo)) {
      data.frame(pos = grid_pos,
                 methylation = approx(.x$pos, .x$mean_meth, xout = grid_pos, rule = 2)$y)
    } else {
      data.frame(pos = grid_pos, methylation = predict(lo, newdata = data.frame(pos = grid_pos)))
    }
  }) %>%
  ungroup()

fig4_trend <- ggplot(smooth_lines, aes(x = pos, y = methylation, color = sample_type)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("WT" = "#2166AC", "KO" = "#D6604D"), name = "") +
  coord_cartesian(xlim = plot_range, ylim = c(0, 100)) +
  theme_pub +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), legend.position = "bottom") +
  labs(y = "Mean methylation (%)")

figure_4 <- grid.arrange(fig4_trend, fig4_tracks, gene_track,
                         ncol = 1, heights = c(0.3, 1.5, 0.3))

ggsave(file.path(OUT_TOX, "Fig4_Tox_CpG_tracks.pdf"), figure_4, width = 12, height = 10, dpi = 300)
ggsave(file.path(OUT_TOX, "Fig4_Tox_CpG_tracks.png"), figure_4, width = 12, height = 10, dpi = 300)
cat("Saved: Fig4_Tox_CpG_tracks\n")

# ==============================================================================
# STEP 12: SUPPLEMENTARY FIGURES
# ==============================================================================

cat("\n======  STEP 12: SUPPLEMENTARY FIGURES  ======\n")

# Load all CpG diff for supplementary
diff_meth_cpg_all <- readRDS(RDS_TOX_DIFF_C)
diff_cpg_df       <- getData(diff_meth_cpg_all) %>% as.data.frame()

# Compute per-CpG methylation from united object
all_cpg_df <- getData(tox_cpg_united) %>%
  as.data.frame() %>%
  mutate(
    WT_mean  = (numCs1/(numCs1+numTs1) + numCs2/(numCs2+numTs2) + numCs3/(numCs3+numTs3)) / 3 * 100,
    KO_mean  = (numCs4/(numCs4+numTs4) + numCs5/(numCs5+numTs5) + numCs6/(numCs6+numTs6)) / 3 * 100
  ) %>%
  left_join(diff_cpg_df %>% dplyr::select(start, meth.diff, qvalue), by = "start") %>%
  mutate(
    significant = start %in% sig_positions,
    meth_diff   = ifelse(is.na(meth.diff), KO_mean - WT_mean, meth.diff),
    qvalue      = replace_na(qvalue, 1)
  )

# --- Supp Fig 1A: Volcano plot ---
supp1a <- ggplot(all_cpg_df, aes(x = meth_diff, y = -log10(qvalue))) +
  geom_point(aes(color = significant), alpha = 0.6, size = 2.5) +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed",
             color = "gray30", linewidth = 0.5) +
  geom_vline(xintercept = c(-20, 20), linetype = "dashed",
             color = "gray30", linewidth = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "#D16666"),
                     labels = c("Non-significant", "Significant DMC"), name = "") +
  annotate("text", x = max(all_cpg_df$meth_diff, na.rm=TRUE) * 0.6,
           y = max(-log10(all_cpg_df$qvalue), na.rm=TRUE) * 0.88,
           label = sprintf("n = %d significant DMCs\n(all hypermethylated)", sum(all_cpg_df$significant)),
           hjust = 0, size = 3.5, fontface = "bold") +
  theme_pub + theme(legend.position = c(0.15, 0.85)) +
  labs(x = "Methylation difference (%) (Tet2KO - WT)",
       y = "-log10(q-value)")

# --- Supp Fig 1B: Histogram of significant DMC meth.diff ---
sig_cpg_data <- all_cpg_df %>% filter(significant)

supp1b <- ggplot(sig_cpg_data, aes(x = meth_diff)) +
  geom_histogram(fill = "grey60", color = "black", bins = 20, alpha = 0.85) +
  geom_vline(xintercept = mean(sig_cpg_data$meth_diff, na.rm = TRUE),
             linetype = "dashed", color = "black", linewidth = 1) +
  annotate("text",
           x = mean(sig_cpg_data$meth_diff, na.rm = TRUE) + 3, y = Inf,
           label = sprintf("Mean: %.1f%%", mean(sig_cpg_data$meth_diff, na.rm=TRUE)),
           hjust = 0, vjust = 1.5, size = 3.5, fontface = "bold") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pub +
  labs(x = "Methylation difference (%)", y = "Number of DMCs")

# --- Supp Fig 1C: WT vs KO scatter ---
supp1c <- ggplot(all_cpg_df, aes(x = WT_mean, y = KO_mean)) +
  geom_point(data = all_cpg_df %>% filter(!significant),
             color = "#D3D3D3", alpha = 0.8, size = 2) +
  geom_point(data = all_cpg_df %>% filter(significant),
             color = "#E57373", size = 2.5, alpha = 0.9) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.8) +
  scale_x_continuous(limits = c(0,100), breaks = seq(0,100,25)) +
  scale_y_continuous(limits = c(0,110), breaks = seq(0,100,25)) +
  theme_pub +
  labs(x = "WT methylation (%)", y = "Tet2KO methylation (%)")

# --- Supp Fig 2: DMR characteristics ---
dmr_data_df <- getData(sig_dmr_tox) %>%
  as.data.frame() %>%
  mutate(width = end - start)

supp2a <- ggplot(dmr_data_df, aes(x = meth.diff)) +
  geom_histogram(fill = "grey60", color = "black", bins = 15, alpha = 0.85) +
  geom_vline(xintercept = mean(dmr_data_df$meth.diff),
             linetype = "dashed", color = "black", linewidth = 1) +
  annotate("text", x = mean(dmr_data_df$meth.diff) + 2, y = Inf,
           label = sprintf("Mean: %.2f%%\n(all hypermethylated)", mean(dmr_data_df$meth.diff)),
           hjust = 0, vjust = 1.5, size = 3.5, fontface = "bold") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_pub +
  labs(x = "Mean methylation difference (%)", y = "Number of DMRs")

supp2b <- ggplot(dmr_data_df, aes(x = width, y = meth.diff)) +
  geom_point(color = "#E74C3C", size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, color = "black",
              fill = "gray80", linewidth = 0.8) +
  stat_cor(method = "pearson", label.x.npc = 0.05, label.y.npc = 0.95, size = 3.5) +
  theme_pub +
  labs(x = "DMR width (bp)", y = "Methylation difference (%)")

# --- Supp Fig 3: Regulatory feature analysis ---
feat_order <- c("enhancer","ctcf","intron","exon","promoter","TSS")
feat_order <- intersect(feat_order, unique(biologically_relevant_dmrs_tox$feature_type))

dmr_annotated_filt <- biologically_relevant_dmrs_tox %>%
  filter(feature_type %in% feat_order) %>%
  mutate(feature_type = factor(feature_type, levels = feat_order))

supp3a <- dmr_annotated_filt %>%
  group_by(feature_type) %>%
  summarise(n = length(unique(peak)), .groups = "drop") %>%
  ggplot(aes(x = reorder(feature_type, -n), y = n, fill = feature_type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5, alpha = 0.9) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = feature_colors) +
  theme_pub + theme(legend.position = "none") +
  labs(x = "Feature type", y = "Number of unique DMRs")

supp3b <- ggplot(dmr_annotated_filt,
                  aes(x = feature_type, y = meth.diff, fill = feature_type)) +
  geom_violin(alpha = 0.7, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.2, fill = "white", outlier.shape = 1, linewidth = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = feature_colors) +
  theme_pub + theme(legend.position = "none",
                    axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Feature type", y = "Methylation difference (%)")

# Save supplementary figures
supp_fig1 <- (supp1a + supp1b) / supp1c +
  plot_annotation(title = "Supplementary Figure: Tox DMC Analysis",
                  theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)))

supp_fig2 <- supp2a + supp2b +
  plot_annotation(title = "Supplementary Figure: Tox DMR Characteristics",
                  theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)))

supp_fig3 <- supp3a + supp3b +
  plot_annotation(title = "Supplementary Figure: Tox Regulatory Element Analysis",
                  theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5)))

ggsave(file.path(OUT_TOX, "Supp_Fig1_DMC_analysis.pdf"),       supp_fig1, width=12, height=10, dpi=300)
ggsave(file.path(OUT_TOX, "Supp_Fig1_DMC_analysis.png"),       supp_fig1, width=12, height=10, dpi=300)
ggsave(file.path(OUT_TOX, "Supp_Fig2_DMR_characteristics.pdf"),supp_fig2, width=10, height=5,  dpi=300)
ggsave(file.path(OUT_TOX, "Supp_Fig2_DMR_characteristics.png"),supp_fig2, width=10, height=5,  dpi=300)
ggsave(file.path(OUT_TOX, "Supp_Fig3_Regulatory_elements.pdf"),supp_fig3, width=10, height=5,  dpi=300)
ggsave(file.path(OUT_TOX, "Supp_Fig3_Regulatory_elements.png"),supp_fig3, width=10, height=5,  dpi=300)
cat("Saved: Supp_Fig1-3\n")

# ==============================================================================
# DONE
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("TOX PIPELINE COMPLETE\n")
cat(strrep("=", 70), "\n\n")
cat("Key outputs in:", OUT_TOX, "\n")
cat("  Final_annotated_biologically_relevant_DMRs_Tox.csv\n")
cat("  Final_annotated_biologically_relevant_DMRs_Tox_coord_adjusted.csv\n")
cat("  Tox_feature_summary.csv\n")
cat("  Tox_unique_DMRs.csv\n")
cat("  Fig1_Tox_CpG_diff_track.pdf/.png\n")
cat("  Fig2_Tox_DMR_landscape.pdf/.png\n")
cat("  Fig3_Tox_methylation_heatmap.pdf\n")
cat("  Fig4_Tox_CpG_tracks.pdf/.png\n")
cat("  Supp_Fig1-3.pdf/.png\n")
