# ============================================================
# CIPN Computational Project
# Script 04: Build GSE185084 DRG sample manifest
# ============================================================

library(GEOquery)

# Re-download/load GEO metadata
gse <- getGEO("GSE185084", GSEMatrix = TRUE)

gse185084 <- gse[[1]]

meta <- pData(gse185084)

# Restrict to DRG sequencing entries
drg_meta <- meta[meta$`tissue:ch1` == "DRGs", ]

# Extract biological sample ID from title
drg_meta$sample_id <- sub(":.*", "", drg_meta$title)

# Keep only relevant columns
manifest_entries <- drg_meta[, c(
  "sample_id",
  "geo_accession",
  "title",
  "treatment:ch1",
  "Sex:ch1",
  "strain:ch1"
)]

# Inspect DRG09 discrepancy
manifest_entries[
  manifest_entries$sample_id == "DRG09",
]
