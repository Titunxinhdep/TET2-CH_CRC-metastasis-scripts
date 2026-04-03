# ==============================================================================
# WGBS COMPLETE ANALYSIS PIPELINE
# TET2KO vs WT CD8+ T cells - Whole Genome Bisulfite Sequencing
# ==============================================================================
# PIPELINE OVERVIEW:
#   STEP 1:  Load and QC methylation data
#   STEP 2:  DMR calling (methylKit)
#   STEP 3:  Coordinate shift GTF (NNNNN-trimmed genome adjustment)
#   STEP 4:  Build annotation object (exon, intron, promoter, TSS, enhancer, CTCF)
#   STEP 5:  Add UTR and CDS to annotation object
#   STEP 6:  Annotate DMRs with all features (annotatePeakInBatch + BioMart)
#   STEP 7:  Gene symbol propagation for enhancer/CTCF
#   STEP 8:  Coordinate shift output CSV (for downstream IGV / figure use)
#   STEP 9:  Overview figures (stacked bar, pie chart, Tox analysis)
#   STEP 10: Supplementary figures (volcano, chr distribution, violin, lollipop)
#   STEP 11: Venn diagram (methylation x expression overlap)
# ==============================================================================

# ==============================================================================
# GLOBAL PATHS  -- edit only this section when moving to a new machine
# ==============================================================================

BASE        <- "/Users/sakatama/Dropbox/Personal/WGBS/Bismark_pipeline"
OUT_MKIT    <- file.path(BASE, "methylKit")
OUT_ALL     <- file.path(OUT_MKIT, "manuscript")

# Raw Bismark coverage files
COVERAGE_FILES <- list(
  file.path(BASE, "Coverage/PRJ76-02_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"),
  file.path(BASE, "Coverage/PRJ76-04_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"),
  file.path(BASE, "Coverage/PRJ76-06_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"),
  file.path(BASE, "Coverage/PRJ76-08_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"),
  file.path(BASE, "Coverage/PRJ76-10_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz"),
  file.path(BASE, "Coverage/PRJ76-12_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz")
)
SAMPLE_IDS  <- list("WT_02", "WT_04", "WT_06", "KO_08", "KO_10", "KO_12")
TREATMENT   <- c(0, 0, 0, 1, 1, 1)   # 0 = WT, 1 = KO

# Reference / annotation files
GTF_ORIGINAL  <- file.path(BASE, "gencode.vM35.primary_assembly.annotation.gtf.gz")
GTF_SHIFTED   <- file.path(BASE, "gencode.vM35.shifted.gtf")
SHIFT_FILE    <- file.path(BASE, "coordinate_shifts_report.csv")
ENHANCER_FILE <- file.path(BASE, "encode_cCRE_enhancers_mm39.bed")
CTCF_FILE     <- file.path(BASE, "encode_cCRE_ctcf_mm39.bed")

# RDS checkpoints (written once, reloaded on subsequent runs)
RDS_METH_OBJ    <- file.path(OUT_MKIT, "meth_object_CD8.rds")
RDS_FILTERED    <- file.path(OUT_MKIT, "meth_object_CD8_lo.count_2_hi.perc_99.9.rds")
RDS_NORMALIZED  <- file.path(OUT_MKIT, "meth_normalized_before_unite.rds")
RDS_UNITED      <- file.path(OUT_MKIT, "meth_unite_CD8.rds")
RDS_DIFFMETH    <- file.path(OUT_MKIT, "unfiltered_diffmeth.rds")
RDS_SIGDMR      <- file.path(OUT_MKIT, "significant_DMRs_qvalue_0.01_meth_cutoff_20.rds")
RDS_ANNOT_OBJ   <- file.path(OUT_MKIT, "annotation_object_complete_with_UTR_CDS.rds")

# Final annotated CSV (canonical output used by all figure scripts)
CSV_FINAL       <- file.path(OUT_ALL, "Final_annotated_biologically_relevant_DMRs.csv")
CSV_FINAL_ADJ   <- file.path(OUT_ALL, "Final_annotated_biologically_relevant_DMRs_coord_adjusted.csv")

# Create output directory if needed
dir.create(OUT_ALL, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# LIBRARIES
# ==============================================================================

suppressPackageStartupMessages({
  library(methylKit)
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(txdbmaker)
  library(ChIPpeakAnno)
  library(biomaRt)
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(VennDiagram)
  library(grid)
  library(gridExtra)
  library(readr)
  library(tidyr)
  library(ggrepel)
  library(RColorBrewer)
})

# ==============================================================================
# STEP 1: LOAD AND QC METHYLATION DATA
# ==============================================================================

cat("\n======  STEP 1: LOAD RAW METHYLATION DATA  ======\n")

if (!file.exists(RDS_METH_OBJ)) {

  # Estimate mean coverage per sample to inform mincov
  get_mean_coverage <- function(file) {
    cov_data <- data.table::fread(file, header = FALSE)
    coverage <- cov_data$V5 + cov_data$V6
    return(mean(coverage))
  }
  mean_cov <- sapply(COVERAGE_FILES, get_mean_coverage)
  names(mean_cov) <- basename(as.character(COVERAGE_FILES))
  cat("Mean coverage per sample:\n"); print(round(mean_cov, 2))

  meth_list <- list()
  for (i in seq_along(COVERAGE_FILES)) {
    tmp <- methRead(
      location  = COVERAGE_FILES[i],
      sample.id = SAMPLE_IDS[i],
      assembly  = "mm39",
      treatment = TREATMENT[i],
      context   = "CpG",
      pipeline  = "bismarkCoverage",
      header    = FALSE,
      mincov    = 2
    )
    meth_list[[i]] <- tmp[[1]]
  }
  meth_obj <- new("methylRawList", meth_list)
  meth_obj@treatment <- TREATMENT

  saveRDS(meth_obj, RDS_METH_OBJ)
  cat("Saved:", RDS_METH_OBJ, "\n")

} else {
  meth_obj <- readRDS(RDS_METH_OBJ)
  cat("Loaded:", RDS_METH_OBJ, "\n")
}

# ==============================================================================
# STEP 2: DMR CALLING
# ==============================================================================

cat("\n======  STEP 2: DMR CALLING  ======\n")

if (!file.exists(RDS_SIGDMR)) {

  filtered_obj  <- filterByCoverage(meth_obj, lo.count = 2, hi.perc = 99.9)
  normalized_obj <- normalizeCoverage(filtered_obj)
  meth_united   <- methylKit::unite(normalized_obj, destrand = FALSE)

  # QC plots
  clusterSamples(meth_united, dist = "correlation", method = "ward", plot = TRUE)
  PCASamples(meth_united)

  # Tile-based DMR calling (1 kb windows)
  tiles     <- tileMethylCounts(meth_united, win.size = 1000, step.size = 1000)
  diff_meth <- calculateDiffMeth(tiles)
  sig_dmr   <- getMethylDiff(diff_meth, difference = 20, qvalue = 0.01)

  cat("Total significant DMRs (|Δmeth|≥20%, q<0.01):", length(sig_dmr), "\n")

  saveRDS(filtered_obj,   RDS_FILTERED)
  saveRDS(normalized_obj, RDS_NORMALIZED)
  saveRDS(meth_united,    RDS_UNITED)
  saveRDS(diff_meth,      RDS_DIFFMETH)
  saveRDS(sig_dmr,        RDS_SIGDMR)

  write.table(as.data.frame(sig_dmr),
              file = file.path(OUT_MKIT, "significant_DMRs_qvalue_0.01_meth_cutoff_20.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

} else {
  sig_dmr   <- readRDS(RDS_SIGDMR)
  diff_meth <- readRDS(RDS_DIFFMETH)
  cat("Loaded:", RDS_SIGDMR, "\n")
  cat("Total significant DMRs:", length(sig_dmr), "\n")
}

# ==============================================================================
# STEP 3: COORDINATE-SHIFT GTF
# (Genome was trimmed of leading NNNNN before alignment; shift back to mm39)
# ==============================================================================

cat("\n======  STEP 3: SHIFT GTF COORDINATES  ======\n")

if (!file.exists(GTF_SHIFTED)) {
  gtf        <- import(GTF_ORIGINAL)
  shift_df   <- read.csv(SHIFT_FILE)
  shift_vec  <- setNames(shift_df$shift, shift_df$chr)

  chr_names     <- as.character(seqnames(gtf))
  shift_amounts <- shift_vec[chr_names]
  shift_amounts[is.na(shift_amounts)] <- 0

  start(gtf) <- start(gtf) - shift_amounts
  end(gtf)   <- end(gtf)   - shift_amounts

  export(gtf, GTF_SHIFTED)
  cat("Shifted GTF saved:", GTF_SHIFTED, "\n")
}

gtf_shifted <- import(GTF_SHIFTED)
txdb        <- makeTxDbFromGRanges(gtf_shifted)
genome(txdb) <- "mm39"
cat("TxDb built from shifted GTF.\n")

# ==============================================================================
# STEP 4: BUILD BASE ANNOTATION OBJECT
# (promoter, TSS, exon, intron[exon-overlap-filtered], enhancer, CTCF)
# ==============================================================================

cat("\n======  STEP 4: BUILD ANNOTATION OBJECT  ======\n")

build_base_annotation <- function(txdb, enhancer_path, ctcf_path, shift_path) {

  shift_df  <- read.csv(shift_path)
  shift_vec <- setNames(shift_df$shift, shift_df$chr)

  apply_shift <- function(gr) {
    seqn   <- as.character(seqnames(gr))
    shifts <- shift_vec[seqn]
    shifts[is.na(shifts)] <- 0
    start(gr) <- start(gr) - shifts
    end(gr)   <- end(gr)   - shifts
    gr
  }

  # 1. Exons
  cat("  Extracting exons...\n")
  exons_by_tx <- exonsBy(txdb, by = "tx")
  txs         <- transcripts(txdb)
  names(exons_by_tx) <- txs$tx_name
  exons_gr    <- unlist(exons_by_tx, use.names = TRUE)
  names(exons_gr) <- rep(names(exons_by_tx), elementNROWS(exons_by_tx))

  valid_tx    <- unique(names(exons_gr))
  tx2gene     <- AnnotationDbi::select(txdb, keys = valid_tx,
                                       columns = c("TXNAME", "GENEID"), keytype = "TXNAME")
  tx2gene_map <- setNames(tx2gene$GENEID, tx2gene$TXNAME)
  exons_gr$gene_id <- tx2gene_map[names(exons_gr)]
  exons_gr$feature <- "exon"

  # 2. Promoters (2 kb upstream, 200 bp downstream)
  cat("  Extracting promoters...\n")
  promoters_gr <- promoters(txdb, upstream = 2000, downstream = 200)
  tx_proms     <- names(promoters_gr)
  proms_df     <- AnnotationDbi::select(txdb, keys = tx_proms,
                                        columns = c("TXNAME", "GENEID"), keytype = "TXNAME")
  proms_df     <- proms_df[match(tx_proms, proms_df$TXNAME), ]
  promoters_gr$gene_id <- proms_df$GENEID
  promoters_gr$feature <- "promoter"

  # 3. TSS (±1 bp)
  cat("  Extracting TSS...\n")
  tss_gr   <- promoters(txdb, upstream = 1, downstream = 1)
  tx_tss   <- names(tss_gr)
  tss_df   <- AnnotationDbi::select(txdb, keys = tx_tss,
                                    columns = c("TXNAME", "GENEID"), keytype = "TXNAME")
  tss_df   <- tss_df[match(tx_tss, tss_df$TXNAME), ]
  tss_gr$gene_id <- tss_df$GENEID
  tss_gr$feature <- "TSS"

  # 4. Introns (remove intron-exon overlaps to avoid false annotation)
  cat("  Extracting introns (filtering exon overlaps)...\n")
  introns_by_tx <- intronsByTranscript(txdb, use.names = TRUE)
  introns_gr    <- unlist(introns_by_tx, use.names = TRUE)
  tx_ids_int    <- rep(names(introns_by_tx), elementNROWS(introns_by_tx))
  names(introns_gr)    <- tx_ids_int
  introns_gr$gene_id   <- tx2gene_map[names(introns_gr)]
  introns_gr$feature   <- "intron"

  overlaps <- findOverlaps(introns_gr, exons_gr)
  if (length(overlaps) > 0) {
    introns_gr <- introns_gr[-unique(queryHits(overlaps))]
  }

  introns_df <- data.frame(
    seqnames = as.character(seqnames(introns_gr)),
    start    = start(introns_gr),
    end      = end(introns_gr),
    strand   = as.character(strand(introns_gr)),
    gene_id  = introns_gr$gene_id,
    tx_id    = names(introns_gr),
    stringsAsFactors = FALSE
  ) %>%
    group_by(seqnames, start, end, strand) %>%
    arrange(desc(!is.na(gene_id))) %>%
    dplyr::slice(1) %>%
    ungroup()

  introns_clean <- makeGRangesFromDataFrame(introns_df, keep.extra.columns = FALSE)
  names(introns_clean)    <- introns_df$tx_id
  introns_clean$gene_id   <- introns_df$gene_id
  introns_clean$feature   <- "intron"

  # 5. Enhancers (ENCODE cCRE, coordinate-shifted)
  cat("  Loading enhancers...\n")
  enhancer_gr <- import(enhancer_path)
  genome(enhancer_gr) <- "mm39"
  enhancer_gr <- apply_shift(enhancer_gr)
  enhancer_gr$feature <- "enhancer"
  names(enhancer_gr)  <- enhancer_gr$name

  # 6. CTCF (ENCODE cCRE, coordinate-shifted)
  cat("  Loading CTCF sites...\n")
  ctcf_gr <- import(ctcf_path)
  genome(ctcf_gr) <- "mm39"
  ctcf_gr <- apply_shift(ctcf_gr)
  ctcf_gr$feature <- "ctcf"
  names(ctcf_gr)  <- ctcf_gr$name

  GRangesList(
    promoters = promoters_gr,
    TSS       = tss_gr,
    exons     = exons_gr,
    introns   = introns_clean,
    enhancers = enhancer_gr,
    ctcf      = ctcf_gr
  )
}

# ==============================================================================
# STEP 5: ADD 5'UTR, 3'UTR AND CDS TO ANNOTATION OBJECT
# (GENCODE labels all UTRs as "UTR"; classify by position relative to CDS)
# ==============================================================================

add_utr_cds <- function(base_annotation, gtf_shifted_gr) {

  cat("  Classifying UTRs (5'/3') from GTF...\n")

  utrs <- gtf_shifted_gr[gtf_shifted_gr$type == "UTR"]
  cds  <- gtf_shifted_gr[gtf_shifted_gr$type == "CDS"]

  utrs_by_tx <- split(utrs, utrs$transcript_id)
  cds_by_tx  <- split(cds,  cds$transcript_id)
  common_txs <- intersect(names(utrs_by_tx), names(cds_by_tx))
  cat("  Transcripts with both UTR and CDS:", length(common_txs), "\n")

  utr3_list <- GRangesList()
  utr5_list <- GRangesList()

  for (tx_id in common_txs) {
    tx_utrs  <- utrs_by_tx[[tx_id]]
    tx_cds   <- cds_by_tx[[tx_id]]
    if (length(tx_utrs) == 0 || length(tx_cds) == 0) next

    cds_start  <- min(start(tx_cds))
    cds_end    <- max(end(tx_cds))
    strand_val <- as.character(strand(tx_cds[1]))

    for (i in seq_along(tx_utrs)) {
      utr       <- tx_utrs[i]
      utr_start <- start(utr)
      utr_end   <- end(utr)

      if (strand_val == "+") {
        if (utr_end <= cds_start) {
          utr$utr_type <- "5UTR"
          utr5_list <- c(utr5_list, GRangesList(utr))
        } else if (utr_start >= cds_end) {
          utr$utr_type <- "3UTR"
          utr3_list <- c(utr3_list, GRangesList(utr))
        }
      } else if (strand_val == "-") {
        if (utr_start >= cds_end) {
          utr$utr_type <- "5UTR"
          utr5_list <- c(utr5_list, GRangesList(utr))
        } else if (utr_end <= cds_start) {
          utr$utr_type <- "3UTR"
          utr3_list <- c(utr3_list, GRangesList(utr))
        }
      }
    }
  }

  utr3 <- unlist(utr3_list)
  utr5 <- unlist(utr5_list)

  clean_meta <- function(gr, feat_label) {
    names(gr) <- gr$transcript_id
    mcols(gr) <- DataFrame(
      gene_id      = gr$gene_id,
      feature      = feat_label,
      feature_type = feat_label
    )
    gr
  }

  # CDS
  cds_ranges <- gtf_shifted_gr[gtf_shifted_gr$type == "CDS"]
  names(cds_ranges) <- cds_ranges$transcript_id
  mcols(cds_ranges) <- DataFrame(
    gene_id      = cds_ranges$gene_id,
    feature      = "CDS",
    feature_type = "CDS"
  )

  annot <- base_annotation
  if (length(utr3) > 0) annot[["3UTR"]] <- clean_meta(utr3, "3UTR")
  if (length(utr5) > 0) annot[["5UTR"]] <- clean_meta(utr5, "5UTR")
  annot[["CDS"]] <- cds_ranges

  cat("  Final annotation features:", paste(names(annot), collapse = ", "), "\n")
  annot
}

# Build or reload annotation object
if (!file.exists(RDS_ANNOT_OBJ)) {
  base_annot      <- build_base_annotation(txdb, ENHANCER_FILE, CTCF_FILE, SHIFT_FILE)
  complete_annot  <- add_utr_cds(base_annot, gtf_shifted)
  saveRDS(complete_annot, RDS_ANNOT_OBJ)
  cat("Saved annotation object:", RDS_ANNOT_OBJ, "\n")
} else {
  complete_annot <- readRDS(RDS_ANNOT_OBJ)
  cat("Loaded annotation object:", RDS_ANNOT_OBJ, "\n")
}

# ==============================================================================
# STEP 6: ANNOTATE DMRs WITH annotatePeakInBatch + BioMart GENE SYMBOLS
# ==============================================================================

cat("\n======  STEP 6: ANNOTATE DMRs  ======\n")

if (!file.exists(CSV_FINAL)) {

  # Prepare GRanges from DMRs
  granges_obj <- as(sig_dmr, "GRanges")
  genome(granges_obj) <- "mm39"

  # Keep only seqlevels present in promoters (representative feature)
  common_seqlevels <- intersect(seqlevels(granges_obj),
                                seqlevels(complete_annot$promoters))
  granges_obj <- keepSeqlevels(granges_obj, common_seqlevels, pruning.mode = "coarse")

  # annotatePeakInBatch for each feature
  annotate_feature <- function(peaks, feat_gr, label) {
    res <- annotatePeakInBatch(peaks, AnnotationData = feat_gr)
    mcols(res)$feature_type <- label
    res
  }

  combined_annotation <- list(
    promoters = annotate_feature(granges_obj, complete_annot$promoters, "promoter"),
    exons     = annotate_feature(granges_obj, complete_annot$exons,     "exon"),
    introns   = annotate_feature(granges_obj, complete_annot$introns,   "intron"),
    TSS       = annotate_feature(granges_obj, complete_annot$TSS,       "TSS"),
    enhancers = annotate_feature(granges_obj, complete_annot$enhancers, "enhancer"),
    ctcf      = annotate_feature(granges_obj, complete_annot$ctcf,      "ctcf"),
    UTR5      = annotate_feature(granges_obj, complete_annot[["5UTR"]], "5UTR"),
    UTR3      = annotate_feature(granges_obj, complete_annot[["3UTR"]], "3UTR"),
    CDS       = annotate_feature(granges_obj, complete_annot$CDS,       "CDS")
  )

  combined_df <- do.call(rbind, lapply(combined_annotation, as.data.frame, row.names = NULL))
  cat("Total annotation rows:", nrow(combined_df), "\n")
  cat("Feature type distribution:\n"); print(table(combined_df$feature_type))

  # BioMart: map transcript IDs → gene symbols
  cat("  Querying BioMart for gene symbols...\n")
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

  combined_df <- combined_df %>%
    mutate(transcript_id_clean = sub("\\.\\d+$", "", feature))

  unique_tx <- unique(combined_df$transcript_id_clean)

  tx2gene_full <- getBM(
    attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                   "external_gene_name", "gene_biotype", "description"),
    filters    = "ensembl_transcript_id",
    values     = unique_tx,
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

  combined_df_annotated <- combined_df %>%
    left_join(tx2gene_agg, by = "transcript_id_clean") %>%
    mutate(gene_symbols = ifelse(is.na(gene_symbols) | gene_symbols == "",
                                 feature, gene_symbols))

  # ==============================================================================
  # STEP 7: PROPAGATE GENE SYMBOLS TO ENHANCER / CTCF VIA OVERLAP
  # ==============================================================================

  cat("  Propagating gene symbols to enhancer/CTCF features...\n")

  gr_all <- makeGRangesFromDataFrame(combined_df_annotated,
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
  final_gr    <- c(retained_gr, regulatory_gr)
  final_df    <- as.data.frame(final_gr) %>% arrange(seqnames, start)

  # Filter for biologically meaningful overlap relationships
  meaningful_annotations <- c("inside", "overlapStart", "overlapEnd", "includeFeature")
  biologically_relevant_dmrs <- final_df %>%
    filter(insideFeature %in% meaningful_annotations)

  cat("Final unique DMRs (biologically relevant):",
      length(unique(biologically_relevant_dmrs$peak)), "\n")

  write.csv(final_df,                   file.path(OUT_ALL, "Final_all_DMRs_annotated_full.csv"),             row.names = FALSE)
  write.csv(biologically_relevant_dmrs, CSV_FINAL,                                                            row.names = FALSE)
  cat("Saved:", CSV_FINAL, "\n")

} else {
  biologically_relevant_dmrs <- read.csv(CSV_FINAL)
  cat("Loaded:", CSV_FINAL, "\n")
  cat("Unique biologically relevant DMRs:", length(unique(biologically_relevant_dmrs$peak)), "\n")
}

# ==============================================================================
# STEP 8: COORDINATE SHIFT OUTPUT CSV (shift back to true mm39 coordinates)
# ==============================================================================

cat("\n======  STEP 8: COORDINATE ADJUSTMENT  ======\n")

if (!file.exists(CSV_FINAL_ADJ)) {

  shift_df <- read.csv(SHIFT_FILE)

  data_adjusted <- biologically_relevant_dmrs %>%
    left_join(shift_df, by = c("seqnames" = "chr")) %>%
    mutate(
      start            = start            + shift,
      end              = end              + shift,
      start_position   = start_position   + shift,
      end_position     = end_position     + shift
    ) %>%
    dplyr::select(-shift) %>%
    mutate(meth_status = ifelse(meth.diff > 0, "hyper", "hypo"))

  # Place meth_status immediately after meth.diff
  idx <- which(names(data_adjusted) == "meth.diff")
  data_adjusted <- data_adjusted[, c(1:idx, ncol(data_adjusted),
                                     (idx+1):(ncol(data_adjusted)-1))]

  write.csv(data_adjusted, CSV_FINAL_ADJ, row.names = FALSE)
  cat("Saved coordinate-adjusted CSV:", CSV_FINAL_ADJ, "\n")

} else {
  data_adjusted <- read.csv(CSV_FINAL_ADJ)
  cat("Loaded:", CSV_FINAL_ADJ, "\n")
}

# ==============================================================================
# COLOR SCHEMES (shared across all figure steps)
# ==============================================================================

feature_colors <- c(
  "enhancer" = "#A8D8A8",
  "promoter" = "#B0C4DE",
  "TSS"      = "#F4A0A0",
  "ctcf"     = "#FFD4A3",
  "intron"   = "#D8BFD8",
  "exon"     = "#C7E9B4",
  "CDS"      = "#FAD7A0",
  "3UTR"     = "#F9E79F",
  "5UTR"     = "#A9CCE3"
)

methylation_colors <- c(
  "Hypermethylated" = "#E6A8A8",
  "Hypomethylated"  = "#A8D5D8",
  "Hyper-rDMRs"     = "#E6A8A8",
  "Hypo-rDMRs"      = "#A8D5D8"
)

theme_wgbs <- theme_classic(base_size = 12) +
  theme(
    plot.title    = element_text(size = 14, face = "bold"),
    axis.title    = element_text(size = 12, face = "bold"),
    axis.text     = element_text(size = 11, color = "black"),
    legend.title  = element_text(size = 11, face = "bold"),
    legend.text   = element_text(size = 10),
    panel.border  = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

# ==============================================================================
# HELPER: deduplicate DMRs (one row per peak)
# ==============================================================================

prepare_unique_dmrs <- function(df) {
  df %>%
    mutate(
      methylation_status = ifelse(meth.diff > 0, "Hyper-rDMRs", "Hypo-rDMRs"),
      feature_clean = case_when(
        feature_type == "enhancer" ~ "Enhancer",
        feature_type == "promoter" ~ "Promoter",
        feature_type == "TSS"      ~ "TSS",
        feature_type == "ctcf"     ~ "CTCF",
        feature_type == "intron"   ~ "Intron",
        feature_type == "exon"     ~ "Exon",
        feature_type == "CDS"      ~ "CDS",
        feature_type == "3UTR"     ~ "3'UTR",
        feature_type == "5UTR"     ~ "5'UTR",
        TRUE ~ feature_type
      ),
      region_category = case_when(
        feature_type %in% c("enhancer","promoter","TSS","ctcf") ~ "Regulatory",
        feature_type %in% c("intron","exon","CDS","3UTR","5UTR") ~ "Genic",
        TRUE ~ "Other"
      )
    ) %>%
    arrange(peak, desc(region_category == "Regulatory"), feature_type) %>%
    group_by(peak) %>%
    dplyr::slice(1) %>%
    ungroup()
}

dmr_unique     <- prepare_unique_dmrs(biologically_relevant_dmrs)
regulatory_dmrs <- length(unique(biologically_relevant_dmrs$peak))
hyper_rdmrs    <- sum(dmr_unique$methylation_status == "Hyper-rDMRs")
hypo_rdmrs     <- sum(dmr_unique$methylation_status == "Hypo-rDMRs")

cat("\n=== KEY STATISTICS ===\n")
cat("Total unique regulatory DMRs:", regulatory_dmrs, "\n")
cat("Hyper-rDMRs:", hyper_rdmrs, sprintf("(%.1f%%)\n", hyper_rdmrs/regulatory_dmrs*100))
cat("Hypo-rDMRs:",  hypo_rdmrs,  sprintf("(%.1f%%)\n", hypo_rdmrs/regulatory_dmrs*100))

# ==============================================================================
# STEP 9: OVERVIEW FIGURES
# ==============================================================================

cat("\n======  STEP 9: OVERVIEW FIGURES  ======\n")

# --- Feature distribution matrix (hyper/hypo stacked bar) ---

feature_dist <- biologically_relevant_dmrs %>%
  mutate(methylation_status = ifelse(meth.diff > 0, "Hypermethylated", "Hypomethylated")) %>%
  group_by(feature_type, methylation_status) %>%
  summarise(unique_dmrs = length(unique(peak)), .groups = "drop")

feature_order_plot <- c("enhancer", "promoter", "TSS", "ctcf", "intron", "exon", "CDS", "5UTR", "3UTR")
feature_types      <- intersect(feature_order_plot, unique(feature_dist$feature_type))

meth_matrix <- matrix(0, nrow = 2, ncol = length(feature_types),
                      dimnames = list(c("Hypomethylated","Hypermethylated"), feature_types))
for (i in seq_len(nrow(feature_dist))) {
  r <- feature_dist$methylation_status[i]
  c <- feature_dist$feature_type[i]
  if (c %in% feature_types) meth_matrix[r, c] <- feature_dist$unique_dmrs[i]
}
meth_matrix <- meth_matrix[, colSums(meth_matrix) > 0, drop = FALSE]

pdf(file.path(OUT_ALL, "Fig_Feature_Distribution_Stacked.pdf"), width = 12, height = 8)
par(mar = c(7, 5, 5, 8), bg = "white")
bp <- barplot(meth_matrix,
              beside = FALSE,
              col    = methylation_colors[rownames(meth_matrix)],
              main   = "DMRs by Feature Type and Methylation Status",
              ylab   = "Number of Unique DMRs",
              border = "black", lwd = 1.2, las = 2,
              cex.names = 1.0, cex.lab = 1.2, cex.main = 1.4,
              legend.text = FALSE, axes = TRUE)

total_counts <- colSums(meth_matrix)
text(x = bp, y = total_counts + max(total_counts) * 0.03,
     labels = total_counts, cex = 1.1, xpd = TRUE, font = 2)

for (i in seq_len(ncol(meth_matrix))) {
  hypo_n   <- meth_matrix["Hypomethylated", i]
  hyper_n  <- meth_matrix["Hypermethylated", i]
  tot      <- total_counts[i]
  if (hypo_n  > 0 && (hypo_n/tot*100)  >= 8)
    text(bp[i], hypo_n/2,            paste0(round(hypo_n/tot*100,1),"%"),  cex=0.9, font=2)
  if (hyper_n > 0 && (hyper_n/tot*100) >= 8)
    text(bp[i], hypo_n + hyper_n/2,  paste0(round(hyper_n/tot*100,1),"%"), cex=0.9, font=2)
}

legend("topright",
       legend = c("Hypermethylated","Hypomethylated"),
       fill   = methylation_colors[c("Hypermethylated","Hypomethylated")],
       border = "black", bty = "n", cex = 1.2)
mtext(paste0("Total Regulatory DMRs: ", regulatory_dmrs,
             "  |  Hyper: ", hyper_rdmrs, "  |  Hypo: ", hypo_rdmrs),
      side = 1, line = 5, cex = 1.0)
dev.off()
cat("Saved: Fig_Feature_Distribution_Stacked.pdf\n")

# --- Hyper vs Hypo pie chart ---

hyper_hypo_data <- data.frame(
  Category   = c("Hypermethylated","Hypomethylated"),
  Count      = c(hyper_rdmrs, hypo_rdmrs),
  Percentage = c(round(hyper_rdmrs/regulatory_dmrs*100,1),
                 round(hypo_rdmrs/regulatory_dmrs*100,1))
)

pie_plot <- ggplot(hyper_hypo_data, aes(x = "", y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 1, color = "black", linewidth = 0.5) +
  coord_polar("y", start = 0) +
  geom_text(aes(label = paste0(Category, "\n", Count, " DMRs\n(", Percentage, "%)")),
            position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = methylation_colors) +
  theme_void() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold")) +
  labs(title = "Hypermethylated vs Hypomethylated Regulatory DMRs",
       subtitle = paste0("Total: ", regulatory_dmrs, " unique DMRs"))

ggsave(file.path(OUT_ALL, "Fig_Hyper_vs_Hypo_Pie.pdf"), pie_plot, width=7, height=6, dpi=300)
ggsave(file.path(OUT_ALL, "Fig_Hyper_vs_Hypo_Pie.png"), pie_plot, width=7, height=6, dpi=300)
cat("Saved: Fig_Hyper_vs_Hypo_Pie\n")

# ==============================================================================
# STEP 10: SUPPLEMENTARY FIGURES
# ==============================================================================

cat("\n======  STEP 10: SUPPLEMENTARY FIGURES  ======\n")

all_dmrs     <- readRDS(RDS_DIFFMETH)
all_dmrs_df  <- as.data.frame(getData(all_dmrs))

# --- Volcano plot (all DMRs before filtering) ---

volcano_data <- all_dmrs_df %>%
  mutate(
    neg_log10_qval = -log10(qvalue),
    significance   = case_when(
      meth.diff >=  20 & qvalue < 0.01 ~ "Hypermethylated",
      meth.diff <= -20 & qvalue < 0.01 ~ "Hypomethylated",
      TRUE ~ "Not significant"
    )
  )

volcano_plot <- ggplot(volcano_data, aes(x = meth.diff, y = neg_log10_qval)) +
  geom_hline(yintercept = -log10(0.01), linetype="dashed", color="black", linewidth=0.8) +
  geom_vline(xintercept = c(-20, 20),   linetype="dashed", color="black", linewidth=0.8) +
  geom_point(data = subset(volcano_data, significance == "Not significant"),
             aes(color=significance), size=0.8, alpha=0.3) +
  geom_point(data = subset(volcano_data, significance == "Hypomethylated"),
             aes(color=significance), size=1.2, alpha=0.8) +
  geom_point(data = subset(volcano_data, significance == "Hypermethylated"),
             aes(color=significance), size=1.2, alpha=0.8) +
  scale_color_manual(values = c("Hypermethylated"="#D16666",
                                "Hypomethylated" ="#66B1B9",
                                "Not significant"="grey60"), name="") +
  scale_x_continuous(breaks=seq(-100,100,20), limits=c(-100,100)) +
  guides(color = guide_legend(override.aes = list(size=4, alpha=1))) +
  theme_wgbs +
  theme(legend.position = c(0.15,0.85),
        legend.background = element_rect(fill="white", color="black")) +
  labs(title    = "Volcano plot of all DMRs",
       subtitle = "Thresholds: |Δmeth| ≥ 20%, q < 0.01",
       x = "Methylation difference (%)",
       y = "-log10(q-value)")

ggsave(file.path(OUT_ALL, "Supp_Volcano_Plot.pdf"), volcano_plot, width=8, height=6, dpi=300)
ggsave(file.path(OUT_ALL, "Supp_Volcano_Plot.png"), volcano_plot, width=8, height=6, dpi=300)
cat("Saved: Supp_Volcano_Plot\n")

# --- Chromosome distribution (significant DMRs only) ---

sig_dmr_df <- getData(sig_dmr) %>%
  as.data.frame() %>%
  mutate(
    methylation_status = ifelse(meth.diff > 0, "Hypermethylated", "Hypomethylated"),
    chr = factor(chr, levels = paste0("chr", c(1:19, "X","Y")))
  ) %>%
  filter(!is.na(chr))

chr_data_all <- sig_dmr_df %>%
  group_by(chr, methylation_status) %>%
  summarise(n = n(), .groups="drop")

chr_plot <- ggplot(chr_data_all, aes(x=chr, y=n, fill=methylation_status)) +
  geom_bar(stat="identity", position="stack", color="black", linewidth=0.3) +
  geom_text(data = chr_data_all %>% group_by(chr) %>% summarise(total=sum(n)),
            aes(x=chr, y=total, label=total, fill=NULL),
            vjust=-0.4, size=2.8, fontface="bold") +
  scale_fill_manual(values = c("Hypermethylated"="#E6A8A8","Hypomethylated"="#A8D5D8"), name="") +
  scale_y_continuous(expand = expansion(mult=c(0,0.12))) +
  theme_wgbs +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=9),
        legend.position = "top") +
  labs(title    = "Chromosome distribution of all significant DMRs",
       subtitle = sprintf("Total: %s DMRs", format(sum(chr_data_all$n), big.mark=",")),
       x = "Chromosome", y = "Number of DMRs")

ggsave(file.path(OUT_ALL, "Supp_Chromosome_Distribution.pdf"), chr_plot, width=10, height=6, dpi=300)
ggsave(file.path(OUT_ALL, "Supp_Chromosome_Distribution.png"), chr_plot, width=10, height=6, dpi=300)
cat("Saved: Supp_Chromosome_Distribution\n")

# --- Methylation change distribution by feature (violin) ---

meth_dist_plot <- dmr_unique %>%
  filter(feature_clean %in% c("Enhancer","Promoter","TSS","CTCF","Intron","Exon")) %>%
  mutate(feature_clean = factor(feature_clean,
                                levels=c("Enhancer","Promoter","TSS","CTCF","Intron","Exon"))) %>%
  ggplot(aes(x=feature_clean, y=meth.diff, fill=feature_clean)) +
  geom_hline(yintercept=c(-20,20), linetype="dashed", color="grey50", linewidth=0.5) +
  geom_hline(yintercept=0,         linetype="solid",  color="black",  linewidth=0.3) +
  geom_violin(alpha=0.7, color="black", linewidth=0.3, scale="width") +
  geom_boxplot(width=0.2, fill="white", outlier.shape=NA, color="black", linewidth=0.3) +
  stat_summary(fun=median, geom="point", shape=23, size=3, fill="red", color="black") +
  scale_fill_brewer(palette="Set2") +
  scale_y_continuous(breaks=seq(-100,100,20)) +
  theme_wgbs +
  theme(axis.text.x=element_text(angle=45,hjust=1), legend.position="none") +
  labs(title    = "Distribution of methylation changes across genomic features",
       subtitle = "Red diamond = median; dashed lines = ±20% threshold",
       x = "", y = "Methylation difference (%)")

ggsave(file.path(OUT_ALL, "Supp_Methylation_Distribution_Violin.pdf"), meth_dist_plot, width=8, height=6, dpi=300)
ggsave(file.path(OUT_ALL, "Supp_Methylation_Distribution_Violin.png"), meth_dist_plot, width=8, height=6, dpi=300)
cat("Saved: Supp_Methylation_Distribution_Violin\n")

# --- Top 20 genes by DMR count (lollipop) ---

top_genes <- biologically_relevant_dmrs %>%
  filter(!is.na(gene_symbols), gene_symbols != "") %>%
  group_by(gene_symbols) %>%
  summarise(
    n_dmrs         = length(unique(peak)),
    mean_meth_diff = mean(meth.diff),
    .groups = "drop"
  ) %>%
  filter(n_dmrs >= 3) %>%
  arrange(desc(n_dmrs)) %>%
  head(20) %>%
  mutate(direction = ifelse(mean_meth_diff > 0, "Hypermethylated", "Hypomethylated"),
         gene_symbols = reorder(gene_symbols, n_dmrs))

write.csv(as.data.frame(top_genes), file.path(OUT_ALL, "Top_20_genes_by_DMR_count.csv"), row.names=FALSE)

lollipop_plot <- ggplot(top_genes, aes(x=gene_symbols, y=n_dmrs)) +
  geom_segment(aes(x=gene_symbols, xend=gene_symbols, y=0, yend=n_dmrs),
               color="grey70", linewidth=0.8) +
  geom_point(aes(size=abs(mean_meth_diff), fill=direction),
             shape=21, color="black", stroke=0.5) +
  scale_fill_manual(values=c("Hypermethylated"="#E6A8A8","Hypomethylated"="#A8D5D8"), name="Direction") +
  scale_size_continuous(name="Mean |Δmeth| (%)", range=c(4,12)) +
  coord_flip() +
  theme_wgbs +
  theme(axis.text.y=element_text(size=10, face="italic"),
        legend.position="right") +
  labs(title    = "Top 20 genes with multiple regulatory DMRs",
       subtitle = "Genes with ≥3 DMRs; point size = mean methylation change",
       x = "", y = "Number of DMRs")

ggsave(file.path(OUT_ALL, "Supp_Top_Genes_Lollipop.pdf"), lollipop_plot, width=10, height=8, dpi=300)
ggsave(file.path(OUT_ALL, "Supp_Top_Genes_Lollipop.png"), lollipop_plot, width=10, height=8, dpi=300)
cat("Saved: Supp_Top_Genes_Lollipop\n")

# --- Feature co-occurrence heatmap ---

feature_binary <- biologically_relevant_dmrs %>%
  select(peak, feature_type) %>%
  distinct() %>%
  mutate(present = 1) %>%
  tidyr::pivot_wider(names_from=feature_type, values_from=present, values_fill=0)

cooccur_features <- c("enhancer","promoter","TSS","ctcf","intron","exon")
cooccur_features <- intersect(cooccur_features, names(feature_binary))
cooccur_mat      <- matrix(0, nrow=length(cooccur_features), ncol=length(cooccur_features),
                           dimnames=list(cooccur_features, cooccur_features))
for (i in seq_along(cooccur_features))
  for (j in seq_along(cooccur_features))
    cooccur_mat[i,j] <- sum(feature_binary[[cooccur_features[i]]] == 1 &
                            feature_binary[[cooccur_features[j]]] == 1)

cooccur_df <- as.data.frame(cooccur_mat) %>%
  tibble::rownames_to_column("Feature1") %>%
  tidyr::pivot_longer(-Feature1, names_to="Feature2", values_to="count") %>%
  mutate(Feature1 = factor(Feature1, levels=cooccur_features),
         Feature2 = factor(Feature2, levels=rev(cooccur_features)))

cooccur_plot <- ggplot(cooccur_df, aes(x=Feature1, y=Feature2, fill=count)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=count), color="white", fontface="bold", size=4) +
  scale_fill_gradient(low="#3498DB", high="#E74C3C", name="DMRs") +
  theme_minimal(base_size=12) +
  theme(axis.text.x=element_text(angle=45,hjust=1), panel.grid=element_blank(),
        plot.title=element_text(face="bold")) +
  labs(title="Feature co-occurrence matrix",
       subtitle="Number of DMRs overlapping each pair of features",
       x="", y="")

ggsave(file.path(OUT_ALL, "Supp_Feature_Cooccurrence.pdf"), cooccur_plot, width=7, height=6, dpi=300)
ggsave(file.path(OUT_ALL, "Supp_Feature_Cooccurrence.png"), cooccur_plot, width=7, height=6, dpi=300)
cat("Saved: Supp_Feature_Cooccurrence\n")

# ==============================================================================
# STEP 11: VENN DIAGRAM (methylation × expression overlap)
# ==============================================================================
# Requires the manually curated CSV with 4 gene-list columns:
#   Upregulated_genes_Vav1Cre_vs_WT
#   Hypomethylated_Vav1Cre_vs_WT
#   Downregulated_genes_Vav1Cre_vs_WT
#   Hypermethylated_Vav1Cre_vs_WT
# ==============================================================================

cat("\n======  STEP 11: VENN DIAGRAM  ======\n")

VENN_INPUT <- file.path(OUT_ALL, "gene_overlap_methylation_expression.csv")
VENN_OUTPUT <- file.path(OUT_ALL, "Venn_methylation_expression_overlap.tiff")

if (file.exists(VENN_INPUT)) {

  venn_data <- read_csv(VENN_INPUT, show_col_types = FALSE)

  extract_genes <- function(col) {
    x <- col[!is.na(col) & col != ""]
    unique(x)
  }

  upregulated      <- extract_genes(venn_data$Upregulated_genes_Vav1Cre_vs_WT)
  hypomethylated   <- extract_genes(venn_data$Hypomethylated_Vav1Cre_vs_WT)
  downregulated    <- extract_genes(venn_data$Downregulated_genes_Vav1Cre_vs_WT)
  hypermethylated  <- extract_genes(venn_data$Hypermethylated_Vav1Cre_vs_WT)

  cat("Gene set sizes for Venn:\n")
  cat("  Upregulated:    ", length(upregulated),    "\n")
  cat("  Hypomethylated: ", length(hypomethylated),  "\n")
  cat("  Downregulated:  ", length(downregulated),   "\n")
  cat("  Hypermethylated:", length(hypermethylated),  "\n")

  make_venn <- function(set1, set2, label1, label2) {
    venn.diagram(
      x        = setNames(list(set1, set2), c(label1, label2)),
      filename = NULL,
      height = 3000, width = 3000, resolution = 300,
      lwd = 2, col = c("black","black"),
      fill = c("lightblue","lightcoral"),
      cex = 0, cat.cex = 0, scaled = TRUE, alpha = 0.7
    )
  }

  tiff(VENN_OUTPUT, width=6000, height=3000, res=300)
  grid.arrange(
    grobTree(make_venn(upregulated,   hypomethylated,  "Upregulated",  "Hypomethylated")),
    grobTree(make_venn(downregulated, hypermethylated, "Downregulated","Hypermethylated")),
    ncol=2
  )
  dev.off()
  cat("Saved:", VENN_OUTPUT, "\n")

} else {
  cat("Venn input file not found:", VENN_INPUT, "\n")
  cat("Place the gene-overlap CSV at that path and rerun STEP 11.\n")
}

# ==============================================================================
# DONE
# ==============================================================================

cat("\n", strrep("=", 70), "\n")
cat("PIPELINE COMPLETE\n")
cat(strrep("=", 70), "\n\n")
cat("Key output files in:", OUT_ALL, "\n")
cat("  Final_annotated_biologically_relevant_DMRs.csv\n")
cat("  Final_annotated_biologically_relevant_DMRs_coord_adjusted.csv\n")
cat("  Fig_Feature_Distribution_Stacked.pdf\n")
cat("  Fig_Hyper_vs_Hypo_Pie.pdf/.png\n")
cat("  Supp_Volcano_Plot.pdf/.png\n")
cat("  Supp_Chromosome_Distribution.pdf/.png\n")
cat("  Supp_Methylation_Distribution_Violin.pdf/.png\n")
cat("  Supp_Top_Genes_Lollipop.pdf/.png\n")
cat("  Supp_Feature_Cooccurrence.pdf/.png\n")
cat("  Venn_methylation_expression_overlap.tiff\n")
