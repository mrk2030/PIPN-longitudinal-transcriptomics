# ============================================================
# CIPN Computational Project
# GSE185084 DRG Pathway Analysis
# Script 07
# ============================================================

# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(fgsea)
library(msigdbr)

# ------------------------------------------------------------
# 2. Read DESeq2 results generated in Script 06
# ------------------------------------------------------------

d4 <- read.csv(
  "results/GSE185084_DE_Paclitaxel_vs_Vehicle_day4.csv",
  row.names = 1
)

d7 <- read.csv(
  "results/GSE185084_DE_Paclitaxel_vs_Vehicle_day7.csv",
  row.names = 1
)

d11 <- read.csv(
  "results/GSE185084_DE_Paclitaxel_vs_Vehicle_day11.csv",
  row.names = 1
)

d23 <- read.csv(
  "results/GSE185084_DE_Paclitaxel_vs_Vehicle_day23.csv",
  row.names = 1
)

# ------------------------------------------------------------
# 3. Build ranked gene vectors for GSEA
# ------------------------------------------------------------

make_rank <- function(x) {
  
  z <- x$stat
  names(z) <- rownames(x)
  
  # Remove missing statistics
  z <- z[!is.na(z)]
  
  # Remove duplicate gene names, if any
  z <- z[!duplicated(names(z))]
  
  # Rank from strongest positive to strongest negative statistic
  sort(z, decreasing = TRUE)
}

rank_day4  <- make_rank(d4)
rank_day7  <- make_rank(d7)
rank_day11 <- make_rank(d11)
rank_day23 <- make_rank(d23)

# Check number of ranked genes
c(
  day4  = length(rank_day4),
  day7  = length(rank_day7),
  day11 = length(rank_day11),
  day23 = length(rank_day23)
)

# ------------------------------------------------------------
# 4. Load mouse Hallmark pathways
# ------------------------------------------------------------

hallmark_mouse <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)

hallmark_pathways <- split(
  hallmark_mouse$gene_symbol,
  hallmark_mouse$gs_name
)

# Check number of Hallmark pathways
length(hallmark_pathways)

# ------------------------------------------------------------
# 5. Define reusable GSEA function
# ------------------------------------------------------------

run_fgsea <- function(rank_vector, pathways) {
  
  x <- fgsea(
    pathways = pathways,
    stats = rank_vector,
    minSize = 10,
    maxSize = 500
  )
  
  x <- x[order(x$padj), ]
  
  return(x)
}

# ------------------------------------------------------------
# 6. Run Hallmark GSEA at all four time points
# ------------------------------------------------------------

set.seed(123)

fgsea_day4_hallmark <- run_fgsea(
  rank_day4,
  hallmark_pathways
)

fgsea_day7_hallmark <- run_fgsea(
  rank_day7,
  hallmark_pathways
)

fgsea_day11_hallmark <- run_fgsea(
  rank_day11,
  hallmark_pathways
)

fgsea_day23_hallmark <- run_fgsea(
  rank_day23,
  hallmark_pathways
)

# ------------------------------------------------------------
# 7. Count significant Hallmark pathways
# ------------------------------------------------------------

significant_pathway_counts <- c(
  day4 = sum(
    fgsea_day4_hallmark$padj < 0.05,
    na.rm = TRUE
  ),
  
  day7 = sum(
    fgsea_day7_hallmark$padj < 0.05,
    na.rm = TRUE
  ),
  
  day11 = sum(
    fgsea_day11_hallmark$padj < 0.05,
    na.rm = TRUE
  ),
  
  day23 = sum(
    fgsea_day23_hallmark$padj < 0.05,
    na.rm = TRUE
  )
)

significant_pathway_counts

# ------------------------------------------------------------
# 8. Save individual GSEA result tables
# ------------------------------------------------------------

save_fgsea_csv <- function(x, filename) {
  
  df <- as.data.frame(x)
  
  # Convert leadingEdge list-column to text
  df$leadingEdge <- vapply(
    df$leadingEdge,
    function(z) paste(z, collapse = ";"),
    character(1)
  )
  
  write.csv(
    df,
    filename,
    row.names = FALSE
  )
}

save_fgsea_csv(
  fgsea_day4_hallmark,
  "results/GSE185084_GSEA_Hallmark_day4.csv"
)

save_fgsea_csv(
  fgsea_day7_hallmark,
  "results/GSE185084_GSEA_Hallmark_day7.csv"
)

save_fgsea_csv(
  fgsea_day11_hallmark,
  "results/GSE185084_GSEA_Hallmark_day11.csv"
)

save_fgsea_csv(
  fgsea_day23_hallmark,
  "results/GSE185084_GSEA_Hallmark_day23.csv"
)

# ------------------------------------------------------------
# 9. Build cross-time Hallmark pathway table
# ------------------------------------------------------------

extract_fgsea <- function(x, day_label) {
  
  df <- as.data.frame(x)[, c(
    "pathway",
    "NES",
    "pval",
    "padj",
    "size"
  )]
  
  colnames(df) <- c(
    "pathway",
    paste0("NES_day", day_label),
    paste0("pval_day", day_label),
    paste0("padj_day", day_label),
    paste0("size_day", day_label)
  )
  
  df
}

g4 <- extract_fgsea(
  fgsea_day4_hallmark,
  4
)

g7 <- extract_fgsea(
  fgsea_day7_hallmark,
  7
)

g11 <- extract_fgsea(
  fgsea_day11_hallmark,
  11
)

g23 <- extract_fgsea(
  fgsea_day23_hallmark,
  23
)

hallmark_timecourse <- Reduce(
  function(x, y) merge(
    x,
    y,
    by = "pathway",
    all = TRUE
  ),
  list(
    g4,
    g7,
    g11,
    g23
  )
)

dim(hallmark_timecourse)

# ------------------------------------------------------------
# 10. Mark pathway significance at each time point
# ------------------------------------------------------------

hallmark_timecourse$sig_day4 <-
  hallmark_timecourse$padj_day4 < 0.05

hallmark_timecourse$sig_day7 <-
  hallmark_timecourse$padj_day7 < 0.05

hallmark_timecourse$sig_day11 <-
  hallmark_timecourse$padj_day11 < 0.05

hallmark_timecourse$sig_day23 <-
  hallmark_timecourse$padj_day23 < 0.05

# ------------------------------------------------------------
# 11. Identify day-11 pathways that resolve or persist
# ------------------------------------------------------------

day11_resolved_pathways <- hallmark_timecourse[
  hallmark_timecourse$sig_day11 &
    !hallmark_timecourse$sig_day23,
]

day11_persistent_pathways <- hallmark_timecourse[
  hallmark_timecourse$sig_day11 &
    hallmark_timecourse$sig_day23,
]

# Sort by statistical significance
day11_resolved_pathways <- day11_resolved_pathways[
  order(day11_resolved_pathways$padj_day11),
]

day11_persistent_pathways <- day11_persistent_pathways[
  order(day11_persistent_pathways$padj_day23),
]

# ------------------------------------------------------------
# 12. Display number of resolving and persistent pathways
# ------------------------------------------------------------

nrow(day11_resolved_pathways)
nrow(day11_persistent_pathways)

# ------------------------------------------------------------
# 13. Inspect key pathway tables
# ------------------------------------------------------------

day11_resolved_pathways[
  ,
  c(
    "pathway",
    "NES_day11",
    "padj_day11",
    "NES_day23",
    "padj_day23"
  )
]

day11_persistent_pathways[
  ,
  c(
    "pathway",
    "NES_day11",
    "padj_day11",
    "NES_day23",
    "padj_day23"
  )
]

# ------------------------------------------------------------
# 14. Save cross-time pathway tables
# ------------------------------------------------------------

write.csv(
  hallmark_timecourse,
  "results/GSE185084_Hallmark_timecourse.csv",
  row.names = FALSE
)

write.csv(
  day11_resolved_pathways,
  "results/GSE185084_Hallmark_day11_resolved.csv",
  row.names = FALSE
)

write.csv(
  day11_persistent_pathways,
  "results/GSE185084_Hallmark_day11_persistent.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Final completion message
# ------------------------------------------------------------

cat(
  "\nScript 07 Hallmark pathway analysis completed successfully.\n"
)# saved 2026-10-04