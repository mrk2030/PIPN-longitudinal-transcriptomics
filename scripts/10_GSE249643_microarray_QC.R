# ============================================================
# CIPN Computational Project
# Script 10: GSE249643 raw-array quality-control screening
# Purpose: determine whether deposited paclitaxel DRG CEL files
# represent independent biological arrays.
# ============================================================

library(GEOquery)
library(oligo)

# Create required folders
dir.create("data_raw/GSE249643_CEL", recursive = TRUE, showWarnings = FALSE)
dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)

# Six DRG samples examined in the candidate external-evaluation dataset
gsm_ids <- c(
  "GSM7956751", "GSM7956752", "GSM7956753",
  "GSM7956754", "GSM7956755", "GSM7956756"
)

# Download sample-level supplementary CEL files
for (gsm in gsm_ids) {
  getGEOSuppFiles(
    gsm,
    makeDirectory = FALSE,
    baseDir = "data_raw/GSE249643_CEL"
  )
}

cel_files <- list.files(
  "data_raw/GSE249643_CEL",
  pattern = "\\.CEL\\.gz$",
  full.names = TRUE
)

if (length(cel_files) != 6) {
  stop("Expected six DRG CEL files for GSE249643.")
}

# Read raw CEL intensities without normalization
raw249643 <- read.celfiles(cel_files)
raw_expr249643 <- exprs(raw249643)

# The first three selected files correspond to the paclitaxel DRG arrays
# in the deposited sample ordering used during the original analysis.
pac_idx <- 1:3

pair_index <- combn(pac_idx, 2)

duplicate_qc <- data.frame(
  comparison = apply(pair_index, 2, function(z) {
    paste0("Paclitaxel_", z[1], "_vs_", z[2])
  }),
  raw_max_absolute_difference = apply(pair_index, 2, function(z) {
    max(abs(raw_expr249643[, z[1]] - raw_expr249643[, z[2]]))
  }),
  raw_correlation = apply(pair_index, 2, function(z) {
    cor(raw_expr249643[, z[1]], raw_expr249643[, z[2]])
  })
)

write.csv(
  duplicate_qc,
  "results/validation/GSE249643_raw_CEL_duplicate_QC.csv",
  row.names = FALSE
)

print(duplicate_qc)

if (all(duplicate_qc$raw_max_absolute_difference == 0)) {
  message(
    "GSE249643 excluded from formal external evaluation: ",
    "the three paclitaxel DRG CEL files are exactly identical ",
    "across all raw probe intensities."
  )
} else {
  warning("Paclitaxel CEL files were not all identical; re-check sample ordering and metadata.")
}
