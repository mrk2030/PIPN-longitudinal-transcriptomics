# ============================================================
# CIPN Computational Project
# Script 02: Download and inspect GSE185084 metadata
# ============================================================

library(GEOquery)

# Download GEO Series metadata
gse <- getGEO("GSE185084", GSEMatrix = TRUE)

# Check how many ExpressionSet objects were returned
length(gse)

# Select the first ExpressionSet
gse185084 <- gse[[1]]

# Extract sample metadata
meta <- pData(gse185084)

# Inspect dimensions
dim(meta)

# View metadata column names
colnames(meta)

# View first few samples
head(meta)

# Save complete metadata table
    write.csv(
      meta,
      file = "metadata/GSE185084_full_metadata.csv",
      row.names = TRUE
    )
    file.exists("metadata/GSE185084_full_metadata.csv")
    unique(meta$`tissue:ch1`)
    unique(meta$`treatment:ch1`)
    unique(meta$`sex:ch1`)
    unique(meta$`strain:ch1`)
    grep("day|time", colnames(meta), value = TRUE, ignore.case = TRUE)
    table(meta$`tissue:ch1`)
    table(meta$`treatment:ch1`)
    table(
      meta$`tissue:ch1`,
      meta$`treatment:ch1`
    )
    grep(
      "sex",
      colnames(meta),
      value = TRUE,
      ignore.case = TRUE
    )
    grep(
      "strain|tissue|treatment",
      colnames(meta),
      value = TRUE,
      ignore.case = TRUE
    )
    drg_meta <- meta[meta$`tissue:ch1` == "DRGs", ]
    dim(drg_meta)
    table(drg_meta$`treatment:ch1`)
    head(
      drg_meta[, c(
        "geo_accession",
        "title",
        "source_name_ch1",
        "treatment:ch1"
      )],
      20
    )
    drg_meta$bio_id <- sub(":.*", "", drg_meta$title)
    head(
      drg_meta[, c(
        "geo_accession",
        "title",
        "bio_id",
        "treatment:ch1"
      )],
      20
    )
    length(unique(drg_meta$bio_id))
    tapply(
      drg_meta$bio_id,
      drg_meta$`treatment:ch1`,
      function(x) length(unique(x))
    )
    table(
      meta[["sex:ch1"]],
      useNA = "ifany"
    )
    colnames(meta)
    table(
      meta[["Sex:ch1"]],
      useNA = "ifany"
    )
    table(
      drg_meta[["Sex:ch1"]],
      useNA = "ifany"
    )
    head(
      drg_meta[, c(
        "geo_accession",
        "title",
        "source_name_ch1",
        "treatment:ch1"
      )],
      20
    )
    drg_meta$bio_id <- sub(":.*", "", drg_meta$title)
    head(
      drg_meta[, c(
        "geo_accession",
        "title",
        "bio_id",
        "treatment:ch1"
      )],
      20
    )
    length(unique(drg_meta$bio_id))
    tapply(
      drg_meta$bio_id,
      drg_meta$`treatment:ch1`,
      function(x) length(unique(x))
    )
    bio_treatment <- unique(
      drg_meta[, c("bio_id", "treatment:ch1")]
    )
    
    table(bio_treatment$bio_id)
    duplicated_bio <- names(
      which(table(bio_treatment$bio_id) > 1)
    )
    
    duplicated_bio
    drg_meta[
      drg_meta$bio_id %in% duplicated_bio,
      c(
        "geo_accession",
        "title",
        "bio_id",
        "treatment:ch1"
      )
    ]
    table(
      drg_meta$bio_id,
      drg_meta$`treatment:ch1`
    )
    drg_meta[
      drg_meta$bio_id == "DRG09",
      c(
        "geo_accession",
        "title",
        "source_name_ch1",
        "treatment:ch1"
      )
    ]
    drg_meta[
      drg_meta$bio_id == "DRG09",
      c(
        "geo_accession",
        "title",
        "source_name_ch1",
        "characteristics_ch1",
        "characteristics_ch1.1",
        "characteristics_ch1.2",
        "characteristics_ch1.3",
        "treatment:ch1"
      )
    ]
    