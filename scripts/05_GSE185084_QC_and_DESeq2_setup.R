# ============================================================
# CIPN Computational Project
# Script 05: GSE185084 DRG QC and DESeq2 setup
# ============================================================

# ------------------------------------------------------------
# 1. Load packages
# ------------------------------------------------------------

library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# ------------------------------------------------------------
# 2. Define project paths
# ------------------------------------------------------------

counts_file <- "data_processed/GSE185084_DRG_raw_counts_clean.csv"

manifest_candidates <- c(
  "metadata/GSE185084_DRG_manifest.csv",
  "data_processed/GSE185084_DRG_manifest.csv",
  "metadata/GSE185084_DRG_sample_manifest.csv",
  "data_processed/GSE185084_DRG_sample_manifest.csv"
)

manifest_file <- manifest_candidates[file.exists(manifest_candidates)][1]

if (is.na(manifest_file)) {
  stop(
    "DRG manifest file not found. Check the output filename from Script 04."
  )
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 3. Load cleaned raw counts
# ------------------------------------------------------------

counts_df <- read.csv(
  counts_file,
  row.names = 1,
  check.names = FALSE
)

counts_mat <- as.matrix(counts_df)

storage.mode(counts_mat) <- "numeric"

if (any(is.na(counts_mat))) {
  stop("Count matrix contains missing values.")
}

if (any(counts_mat < 0)) {
  stop("Count matrix contains negative values.")
}

if (any(counts_mat %% 1 != 0)) {
  stop("Count matrix contains non-integer values.")
}

storage.mode(counts_mat) <- "integer"

cat("Raw cleaned count matrix:\n")
cat(nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n\n")

# ------------------------------------------------------------
# 4. Load and standardize DRG sample manifest
# ------------------------------------------------------------

manifest <- read.csv(
  manifest_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Using manifest:\n", manifest_file, "\n\n")

# Identify sample ID column
sample_col_candidates <- c(
  "sample_id",
  "sample",
  "Sample",
  "sampleID",
  "sample_ID"
)

sample_col <- sample_col_candidates[
  sample_col_candidates %in% colnames(manifest)
][1]

if (is.na(sample_col)) {
  stop("Could not identify sample-ID column in manifest.")
}

rownames(manifest) <- manifest[[sample_col]]

# Identify treatment column
treatment_candidates <- c(
  "treatment",
  "Treatment",
  "treatment_group",
  "group"
)

treatment_col <- treatment_candidates[
  treatment_candidates %in% colnames(manifest)
][1]

if (is.na(treatment_col)) {
  stop("Could not identify treatment column in manifest.")
}

# Identify day column
day_candidates <- c(
  "day",
  "Day",
  "collection_day",
  "timepoint",
  "time_point"
)

day_col <- day_candidates[
  day_candidates %in% colnames(manifest)
][1]

if (is.na(day_col)) {
  stop("Could not identify day column in manifest.")
}

# ------------------------------------------------------------
# 5. Match metadata to count-matrix columns
# ------------------------------------------------------------

missing_meta <- setdiff(colnames(counts_mat), rownames(manifest))

if (length(missing_meta) > 0) {
  stop(
    "Samples present in count matrix but absent from manifest: ",
    paste(missing_meta, collapse = ", ")
  )
}

manifest <- manifest[colnames(counts_mat), , drop = FALSE]

if (!identical(rownames(manifest), colnames(counts_mat))) {
  stop("Metadata and count matrix sample order do not match.")
}

# ------------------------------------------------------------
# 6. Standardize experimental factors
# ------------------------------------------------------------

treatment_raw <- as.character(manifest[[treatment_col]])

# Standardize common treatment labels
treatment_raw[
  grepl("vehicle|control", treatment_raw, ignore.case = TRUE)
] <- "Vehicle"

treatment_raw[
  grepl("paclitaxel|taxol", treatment_raw, ignore.case = TRUE)
] <- "Paclitaxel"

manifest$treatment <- factor(
  treatment_raw,
  levels = c("Vehicle", "Paclitaxel")
)

day_raw <- gsub("[^0-9]", "", as.character(manifest[[day_col]]))

manifest$day <- factor(
  day_raw,
  levels = c("4", "7", "11", "23")
)

if (any(is.na(manifest$treatment))) {
  stop("Unrecognized treatment labels detected.")
}

if (any(is.na(manifest$day))) {
  stop("Unrecognized collection-day labels detected.")
}

cat("Sample distribution:\n")
print(table(manifest$day, manifest$treatment))
cat("\n")

# ------------------------------------------------------------
# 7. Apply exact low-count filter
#    Retain genes with >=10 raw counts in >=4 samples
# ------------------------------------------------------------

keep <- rowSums(counts_mat >= 10) >= 4

counts_filtered <- counts_mat[keep, , drop = FALSE]

cat("Genes before filtering:", nrow(counts_mat), "\n")
cat("Genes after filtering :", nrow(counts_filtered), "\n\n")

if (nrow(counts_filtered) != 19055) {
  warning(
    "Expected 19,055 genes after filtering, but obtained ",
    nrow(counts_filtered),
    ". Verify input count file and filtering."
  )
}

# ------------------------------------------------------------
# 8. Construct DESeq2 object
# ------------------------------------------------------------

coldata <- data.frame(
  day = manifest$day,
  treatment = manifest$treatment,
  row.names = rownames(manifest)
)

dds <- DESeqDataSetFromMatrix(
  countData = counts_filtered,
  colData = coldata,
  design = ~ day * treatment
)

# Explicitly preserve reference levels
dds$day <- relevel(dds$day, ref = "4")
dds$treatment <- relevel(dds$treatment, ref = "Vehicle")

design(dds) <- ~ day * treatment

# ------------------------------------------------------------
# 9. Fit DESeq2 model
# ------------------------------------------------------------

dds <- DESeq(dds)

cat("DESeq2 coefficient names:\n")
print(resultsNames(dds))
cat("\n")

# Expected:
# Intercept
# day_7_vs_4
# day_11_vs_4
# day_23_vs_4
# treatment_Paclitaxel_vs_Vehicle
# day7.treatmentPaclitaxel
# day11.treatmentPaclitaxel
# day23.treatmentPaclitaxel

saveRDS(
  dds,
  "data_processed/GSE185084_DRG_DESeq2_fitted.rds"
)

# ------------------------------------------------------------
# 10. Variance-stabilizing transformation for QC only
# ------------------------------------------------------------

vsd <- vst(
  dds,
  blind = FALSE
)

saveRDS(
  vsd,
  "data_processed/GSE185084_DRG_VST.rds"
)

# ------------------------------------------------------------
# 11. PCA QC
# ------------------------------------------------------------

pca_data <- plotPCA(
  vsd,
  intgroup = c("day", "treatment"),
  returnData = TRUE
)

percentVar <- round(
  100 * attr(pca_data, "percentVar")
)

pca_data$sample <- rownames(pca_data)

write.csv(
  pca_data,
  "results/GSE185084_DRG_PCA_data.csv",
  row.names = FALSE
)

cat(
  "PCA variance explained:\n",
  "PC1 =", percentVar[1], "%\n",
  "PC2 =", percentVar[2], "%\n\n"
)

# Publication-style PCA
pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    shape = treatment
  )
) +
  geom_point(size = 3.2) +
  geom_text_repel(
    data = subset(pca_data, sample == "DRG09"),
    aes(label = sample),
    size = 4,
    box.padding = 0.5,
    point.padding = 0.4,
    show.legend = FALSE
  ) +
  labs(
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    shape = "Treatment"
  ) +
  theme_classic(base_size = 12)

ggsave(
  "figures/GSE185084_DRG_PCA_PUBLICATION.pdf",
  pca_plot,
  width = 7,
  height = 5
)

# ------------------------------------------------------------
# 12. Sample-distance QC
# ------------------------------------------------------------

sample_dist <- dist(t(assay(vsd)))

sample_dist_matrix <- as.matrix(sample_dist)

write.csv(
  sample_dist_matrix,
  "results/GSE185084_DRG_sample_distance_matrix.csv"
)

annotation_col <- data.frame(
  Day = coldata$day,
  Treatment = coldata$treatment
)

rownames(annotation_col) <- rownames(coldata)

pdf(
  "figures/GSE185084_DRG_sample_distance_heatmap_PUBLICATION.pdf",
  width = 8,
  height = 8
)

pheatmap(
  sample_dist_matrix,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  border_color = NA,
  main = "GSE185084 DRG sample distances"
)

dev.off()

# ------------------------------------------------------------
# 13. Sample-correlation QC
# ------------------------------------------------------------

sample_cor <- cor(
  assay(vsd),
  method = "pearson"
)

write.csv(
  sample_cor,
  "results/GSE185084_DRG_sample_correlations.csv"
)

drg09_cor <- sample_cor[
  "DRG09",
  setdiff(colnames(sample_cor), "DRG09")
]

# ------------------------------------------------------------
# 14. Export QC summary
# ------------------------------------------------------------

qc_summary <- data.frame(
  metric = c(
    "Raw genes",
    "Filtered genes",
    "Samples",
    "PC1 variance (%)",
    "PC2 variance (%)",
    "DRG09 minimum Pearson correlation",
    "DRG09 median Pearson correlation"
  ),
  value = c(
    nrow(counts_mat),
    nrow(counts_filtered),
    ncol(counts_filtered),
    percentVar[1],
    percentVar[2],
    min(drg09_cor, na.rm = TRUE),
    median(drg09_cor, na.rm = TRUE)
  )
)

write.csv(
  qc_summary,
  "results/GSE185084_DRG_QC_summary.csv",
  row.names = FALSE
)

print(qc_summary)

# ------------------------------------------------------------
# 15. Final checks
# ------------------------------------------------------------

stopifnot(
  ncol(dds) == 32,
  nrow(dds) == 19055,
  identical(
    colnames(dds),
    rownames(colData(dds))
  )
)

cat("\n")
cat("Script 05 completed successfully.\n")
cat("All 32 DRG samples retained, including DRG09.\n")
cat("Filtered genes:", nrow(dds), "\n")
cat("DESeq2 fitted object saved successfully.\n")

# ============================================================
# End of Script 05
# ============================================================
