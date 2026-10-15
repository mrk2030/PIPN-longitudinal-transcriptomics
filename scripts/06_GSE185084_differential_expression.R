# =====================================================================
# CIPN Computational Project
# GSE185084
#
# Script 06: Differential expression and longitudinal trajectory analysis
#
# Primary tissue: DRG
# Comparison: Paclitaxel vs Vehicle
# Time points: Days 4, 7, 11 and 23
# =====================================================================


# ---------------------------------------------------------------------
# 1. Load required packages
# ---------------------------------------------------------------------

library(DESeq2)
library(ggplot2)


# ---------------------------------------------------------------------
# 2. Load cleaned count data and finalized sample manifest
# ---------------------------------------------------------------------

counts_df <- read.csv(
  "data_processed/GSE185084_DRG_raw_counts_clean.csv",
  check.names = FALSE
)

manifest <- read.csv(
  "metadata/GSE185084_DRG_sample_manifest.csv",
  stringsAsFactors = FALSE
)


# ---------------------------------------------------------------------
# 3. Construct count matrix
# ---------------------------------------------------------------------

count_matrix <- as.matrix(
  counts_df[, -1]
)

rownames(count_matrix) <- counts_df$gene_symbol


# ---------------------------------------------------------------------
# 4. Match metadata order to count-matrix columns
# ---------------------------------------------------------------------

manifest <- manifest[
  match(colnames(count_matrix), manifest$sample_id),
]

stopifnot(
  identical(
    colnames(count_matrix),
    manifest$sample_id
  )
)


# ---------------------------------------------------------------------
# 5. Define experimental factors
# ---------------------------------------------------------------------

manifest$treatment <- factor(
  manifest$treatment,
  levels = c(
    "Vehicle",
    "Paclitaxel"
  )
)

manifest$day <- factor(
  manifest$day,
  levels = c(
    4,
    7,
    11,
    23
  )
)

rownames(manifest) <- manifest$sample_id


# Final alignment check
stopifnot(
  identical(
    colnames(count_matrix),
    rownames(manifest)
  )
)


# ---------------------------------------------------------------------
# 6. Construct DESeq2 dataset
# ---------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = manifest,
  design = ~ day + treatment + day:treatment
)


# ---------------------------------------------------------------------
# 7. Low-count filtering
#
# Retain genes with >=10 raw counts in at least 4 samples.
# Four samples correspond to the minimum biological group size.
# ---------------------------------------------------------------------

keep <- rowSums(
  counts(dds) >= 10
) >= 4

table(keep)

dds_filtered <- dds[
  keep,
]


# Inspect dimensions
dim(dds)
dim(dds_filtered)


# ---------------------------------------------------------------------
# 8. Fit full DESeq2 interaction model
# ---------------------------------------------------------------------

dds_fit <- DESeq(
  dds_filtered
)

design(dds_fit)

resultsNames(dds_fit)


# Save fitted object
saveRDS(
  dds_fit,
  file =
    "data_processed/GSE185084_DRG_DESeq2_fitted.rds"
)


# =====================================================================
# PART A
# DAY-SPECIFIC PACLITAXEL VS VEHICLE CONTRASTS
# =====================================================================


# ---------------------------------------------------------------------
# 9. Day 4
#
# Day 4 is the reference time point.
# Therefore, the treatment main effect directly represents:
#
# Paclitaxel vs Vehicle at Day 4
# ---------------------------------------------------------------------

res_day4 <- results(
  dds_fit,
  name =
    "treatment_Paclitaxel_vs_Vehicle",
  alpha = 0.05
)

summary(res_day4)


# ---------------------------------------------------------------------
# 10. Day 7
# ---------------------------------------------------------------------

res_day7 <- results(
  dds_fit,
  contrast = list(
    c(
      "treatment_Paclitaxel_vs_Vehicle",
      "day7.treatmentPaclitaxel"
    )
  ),
  alpha = 0.05
)

summary(res_day7)


# ---------------------------------------------------------------------
# 11. Day 11
# ---------------------------------------------------------------------

res_day11 <- results(
  dds_fit,
  contrast = list(
    c(
      "treatment_Paclitaxel_vs_Vehicle",
      "day11.treatmentPaclitaxel"
    )
  ),
  alpha = 0.05
)

summary(res_day11)


# ---------------------------------------------------------------------
# 12. Day 23
# ---------------------------------------------------------------------

res_day23 <- results(
  dds_fit,
  contrast = list(
    c(
      "treatment_Paclitaxel_vs_Vehicle",
      "day23.treatmentPaclitaxel"
    )
  ),
  alpha = 0.05
)

summary(res_day23)


# ---------------------------------------------------------------------
# 13. Count FDR-significant genes
# ---------------------------------------------------------------------

sig_day4 <- sum(
  res_day4$padj < 0.05,
  na.rm = TRUE
)

sig_day7 <- sum(
  res_day7$padj < 0.05,
  na.rm = TRUE
)

sig_day11 <- sum(
  res_day11$padj < 0.05,
  na.rm = TRUE
)

sig_day23 <- sum(
  res_day23$padj < 0.05,
  na.rm = TRUE
)


c(
  day4 = sig_day4,
  day7 = sig_day7,
  day11 = sig_day11,
  day23 = sig_day23
)


# Expected from our current analysis:
#
# Day 4  = 0
# Day 7  = 4
# Day 11 = 62
# Day 23 = 16


# ---------------------------------------------------------------------
# 14. Count upregulated and downregulated genes
# ---------------------------------------------------------------------

count_direction <- function(res) {
  
  c(
    
    up = sum(
      res$padj < 0.05 &
        res$log2FoldChange > 0,
      na.rm = TRUE
    ),
    
    down = sum(
      res$padj < 0.05 &
        res$log2FoldChange < 0,
      na.rm = TRUE
    )
    
  )
}


count_direction(res_day4)
count_direction(res_day7)
count_direction(res_day11)
count_direction(res_day23)


# ---------------------------------------------------------------------
# 15. Build compact DE summary table
# ---------------------------------------------------------------------

de_summary <- data.frame(
  
  day = c(
    4,
    7,
    11,
    23
  ),
  
  significant = c(
    sig_day4,
    sig_day7,
    sig_day11,
    sig_day23
  ),
  
  up = c(
    
    sum(
      res_day4$padj < 0.05 &
        res_day4$log2FoldChange > 0,
      na.rm = TRUE
    ),
    
    sum(
      res_day7$padj < 0.05 &
        res_day7$log2FoldChange > 0,
      na.rm = TRUE
    ),
    
    sum(
      res_day11$padj < 0.05 &
        res_day11$log2FoldChange > 0,
      na.rm = TRUE
    ),
    
    sum(
      res_day23$padj < 0.05 &
        res_day23$log2FoldChange > 0,
      na.rm = TRUE
    )
  ),
  
  down = c(
    
    sum(
      res_day4$padj < 0.05 &
        res_day4$log2FoldChange < 0,
      na.rm = TRUE
    ),
    
    sum(
      res_day7$padj < 0.05 &
        res_day7$log2FoldChange < 0,
      na.rm = TRUE
    ),
    
    sum(
      res_day11$padj < 0.05 &
        res_day11$log2FoldChange < 0,
      na.rm = TRUE
    ),
    
    sum(
      res_day23$padj < 0.05 &
        res_day23$log2FoldChange < 0,
      na.rm = TRUE
    )
  )
)

de_summary


write.csv(
  de_summary,
  "results/GSE185084_DE_summary_by_day.csv",
  row.names = FALSE
)


# ---------------------------------------------------------------------
# 16. Save complete day-specific DESeq2 results
# ---------------------------------------------------------------------

write.csv(
  as.data.frame(res_day4),
  "results/GSE185084_DE_Paclitaxel_vs_Vehicle_day4.csv"
)

write.csv(
  as.data.frame(res_day7),
  "results/GSE185084_DE_Paclitaxel_vs_Vehicle_day7.csv"
)

write.csv(
  as.data.frame(res_day11),
  "results/GSE185084_DE_Paclitaxel_vs_Vehicle_day11.csv"
)

write.csv(
  as.data.frame(res_day23),
  "results/GSE185084_DE_Paclitaxel_vs_Vehicle_day23.csv"
)


# =====================================================================
# PART B
# SIGNIFICANT GENES AT THE PEAK AND LATE TIME POINTS
# =====================================================================


# ---------------------------------------------------------------------
# 17. Day-11 significant genes
# ---------------------------------------------------------------------

day11_sig <- as.data.frame(
  res_day11
)

day11_sig$gene <- rownames(
  day11_sig
)

day11_sig <- day11_sig[
  !is.na(day11_sig$padj) &
    day11_sig$padj < 0.05,
]

day11_sig <- day11_sig[
  order(day11_sig$padj),
]

head(
  day11_sig,
  20
)


write.csv(
  day11_sig,
  "results/GSE185084_day11_significant_genes.csv",
  row.names = FALSE
)


# ---------------------------------------------------------------------
# 18. Day-23 significant genes
# ---------------------------------------------------------------------

day23_sig <- as.data.frame(
  res_day23
)

day23_sig$gene <- rownames(
  day23_sig
)

day23_sig <- day23_sig[
  !is.na(day23_sig$padj) &
    day23_sig$padj < 0.05,
]

day23_sig <- day23_sig[
  order(day23_sig$padj),
]

head(
  day23_sig,
  20
)


write.csv(
  day23_sig,
  "results/GSE185084_day23_significant_genes.csv",
  row.names = FALSE
)


# =====================================================================
# PART C
# GLOBAL TREATMENT x TIME INTERACTION
# =====================================================================


# ---------------------------------------------------------------------
# 19. Likelihood-ratio test
#
# Full model:
# day + treatment + day:treatment
#
# Reduced model:
# day + treatment
#
# This tests whether the paclitaxel effect changes over time.
# ---------------------------------------------------------------------

dds_lrt <- DESeq(
  dds_filtered,
  test = "LRT",
  reduced = ~ day + treatment
)

res_interaction <- results(
  dds_lrt,
  alpha = 0.05
)

summary(
  res_interaction
)


interaction_count <- sum(
  res_interaction$padj < 0.05,
  na.rm = TRUE
)

interaction_count


# Expected from current analysis:
#
# 10 significant treatment x time interaction genes


# ---------------------------------------------------------------------
# 20. Extract significant interaction genes
# ---------------------------------------------------------------------

interaction_sig <- as.data.frame(
  res_interaction
)

interaction_sig$gene <- rownames(
  interaction_sig
)

interaction_sig <- interaction_sig[
  !is.na(interaction_sig$padj) &
    interaction_sig$padj < 0.05,
]

interaction_sig <- interaction_sig[
  order(interaction_sig$padj),
]

interaction_sig


write.csv(
  interaction_sig,
  "results/GSE185084_significant_treatment_time_interaction_genes.csv",
  row.names = FALSE
)


# =====================================================================
# PART D
# BUILD LONGITUDINAL GENE TRAJECTORY TABLE
# =====================================================================


# ---------------------------------------------------------------------
# 21. Function for converting DESeq2 results into standardized tables
# ---------------------------------------------------------------------

make_res_df <- function(
    res,
    day_label
) {
  
  df <- as.data.frame(
    res
  )
  
  df$gene <- rownames(
    df
  )
  
  df <- df[
    ,
    c(
      "gene",
      "baseMean",
      "log2FoldChange",
      "pvalue",
      "padj"
    )
  ]
  
  colnames(df) <- c(
    
    "gene",
    
    paste0(
      "baseMean_day",
      day_label
    ),
    
    paste0(
      "log2FC_day",
      day_label
    ),
    
    paste0(
      "pvalue_day",
      day_label
    ),
    
    paste0(
      "padj_day",
      day_label
    )
  )
  
  df
}


# ---------------------------------------------------------------------
# 22. Generate standardized tables
# ---------------------------------------------------------------------

df4 <- make_res_df(
  res_day4,
  4
)

df7 <- make_res_df(
  res_day7,
  7
)

df11 <- make_res_df(
  res_day11,
  11
)

df23 <- make_res_df(
  res_day23,
  23
)


# ---------------------------------------------------------------------
# 23. Merge all four time points
# ---------------------------------------------------------------------

trajectory <- Reduce(
  
  function(
    x,
    y
  ) {
    
    merge(
      x,
      y,
      by = "gene",
      all = TRUE
    )
    
  },
  
  list(
    df4,
    df7,
    df11,
    df23
  )
)


dim(
  trajectory
)

head(
  trajectory
)


# Expected:
#
# approximately 19,055 genes


# =====================================================================
# PART E
# DEFINE LONGITUDINAL RESPONSE CHARACTERISTICS
# =====================================================================


# ---------------------------------------------------------------------
# 24. Statistical significance indicators
# ---------------------------------------------------------------------

trajectory$day11_sig <-
  !is.na(
    trajectory$padj_day11
  ) &
  trajectory$padj_day11 < 0.05


trajectory$day23_sig <-
  !is.na(
    trajectory$padj_day23
  ) &
  trajectory$padj_day23 < 0.05


# ---------------------------------------------------------------------
# 25. Directional consistency
# ---------------------------------------------------------------------

trajectory$same_direction_11_23 <-
  
  sign(
    trajectory$log2FC_day11
  ) ==
  
  sign(
    trajectory$log2FC_day23
  )


# ---------------------------------------------------------------------
# 26. Relative magnitude remaining at Day 23
#
# |LFC Day23| / |LFC Day11|
# ---------------------------------------------------------------------

trajectory$relative_effect_day23 <-
  
  ifelse(
    
    abs(
      trajectory$log2FC_day11
    ) > 0,
    
    abs(
      trajectory$log2FC_day23
    ) /
      abs(
        trajectory$log2FC_day11
      ),
    
    NA
  )


# ---------------------------------------------------------------------
# 27. Mark significant interaction genes
# ---------------------------------------------------------------------

trajectory$interaction_sig <-
  
  trajectory$gene %in%
  interaction_sig$gene


# =====================================================================
# PART F
# FINAL TRAJECTORY CLASSIFICATION
# =====================================================================


# ---------------------------------------------------------------------
# Definitions
#
# Strongly_resolving:
#   Day-11 significant AND
#   Day-23 absolute effect <50% of Day-11 effect.
#
# Residual:
#   Day-11 significant,
#   same direction at Day 23,
#   >=50% effect retained,
#   but Day-23 FDR >=0.05.
#
# Persistent:
#   Day-11 significant,
#   same direction,
#   >=50% effect retained,
#   and Day-23 FDR <0.05.
#
# The 50% threshold is an operational trajectory criterion and should
# later be tested in sensitivity analyses.
# ---------------------------------------------------------------------

trajectory$final_class <-
  "Not_day11_responsive"


# Strongly resolving
trajectory$final_class[
  
  trajectory$day11_sig &
    
    trajectory$relative_effect_day23 <
    0.5
  
] <- "Strongly_resolving"


# Residual
trajectory$final_class[
  
  trajectory$day11_sig &
    
    trajectory$same_direction_11_23 &
    
    trajectory$relative_effect_day23 >=
    0.5 &
    
    !trajectory$day23_sig
  
] <- "Residual"


# Persistent
trajectory$final_class[
  
  trajectory$day11_sig &
    
    trajectory$same_direction_11_23 &
    
    trajectory$relative_effect_day23 >=
    0.5 &
    
    trajectory$day23_sig
  
] <- "Persistent"


# ---------------------------------------------------------------------
# 28. Identify resolving genes that cross zero slightly
#
# These are retained within Strongly_resolving but separately annotated.
# ---------------------------------------------------------------------

trajectory$resolving_overshoot <-
  
  trajectory$day11_sig &
  
  !trajectory$same_direction_11_23 &
  
  trajectory$relative_effect_day23 <
  0.5


# ---------------------------------------------------------------------
# 29. Inspect final trajectory classes
# ---------------------------------------------------------------------

table(
  trajectory$final_class
)


table(
  trajectory$final_class,
  trajectory$interaction_sig
)


sum(
  trajectory$resolving_overshoot
)


# Expected from current analysis:
#
# Not_day11_responsive   18993
# Persistent                 4
# Residual                   8
# Strongly_resolving        50
#
# resolving_overshoot = 9


# ---------------------------------------------------------------------
# 30. Save final trajectory table
# ---------------------------------------------------------------------

write.csv(
  trajectory,
  "results/GSE185084_longitudinal_gene_trajectories_FINAL.csv",
  row.names = FALSE
)


# =====================================================================
# PART G
# INSPECT SPECIAL TEMPORAL GENE SETS
# =====================================================================


# ---------------------------------------------------------------------
# 31. Resolving genes with small sign-crossing / overshoot
# ---------------------------------------------------------------------

resolving_overshoot_genes <- trajectory[
  trajectory$resolving_overshoot,
  c(
    "gene",
    "log2FC_day4",
    "log2FC_day7",
    "log2FC_day11",
    "log2FC_day23",
    "padj_day11",
    "padj_day23",
    "relative_effect_day23",
    "interaction_sig"
  )
]

resolving_overshoot_genes


write.csv(
  resolving_overshoot_genes,
  "results/GSE185084_resolving_overshoot_genes.csv",
  row.names = FALSE
)


# ---------------------------------------------------------------------
# 32. Interaction genes not significant at Day 11
# ---------------------------------------------------------------------

interaction_non_day11 <- trajectory[
  trajectory$interaction_sig &
    trajectory$final_class ==
    "Not_day11_responsive",
  c(
    "gene",
    "log2FC_day4",
    "log2FC_day7",
    "log2FC_day11",
    "log2FC_day23",
    "padj_day11",
    "padj_day23"
  )
]

interaction_non_day11


write.csv(
  interaction_non_day11,
  "results/GSE185084_interaction_genes_not_day11_significant.csv",
  row.names = FALSE
)


# =====================================================================
# PART H
# LONGITUDINAL TRAJECTORY FIGURE
# =====================================================================


# ---------------------------------------------------------------------
# 33. Select the 62 Day-11 responsive genes
# ---------------------------------------------------------------------

plot_df <- trajectory[
  trajectory$final_class %in%
    c(
      "Strongly_resolving",
      "Residual",
      "Persistent"
    ),
  c(
    "gene",
    "final_class",
    "log2FC_day4",
    "log2FC_day7",
    "log2FC_day11",
    "log2FC_day23"
  )
]


# ---------------------------------------------------------------------
# 34. Convert to long format
# ---------------------------------------------------------------------

plot_long <- reshape(
  
  plot_df,
  
  varying = c(
    "log2FC_day4",
    "log2FC_day7",
    "log2FC_day11",
    "log2FC_day23"
  ),
  
  v.names = "log2FC",
  
  timevar = "day",
  
  times = c(
    4,
    7,
    11,
    23
  ),
  
  direction = "long"
)


plot_long$day <- as.numeric(
  as.character(
    plot_long$day
  )
)


head(
  plot_long
)


# ---------------------------------------------------------------------
# 35. Generate trajectory plot
# ---------------------------------------------------------------------

trajectory_plot <- ggplot(
  
  plot_long,
  
  aes(
    x = day,
    y = log2FC,
    group = gene
  )
  
) +
  
  geom_line(
    alpha = 0.35
  ) +
  
  geom_point(
    size = 1.5
  ) +
  
  facet_wrap(
    ~ final_class,
    scales = "free_y",
    labeller = as_labeller(c(
      "Persistent" = "Persistent",
      "Residual" = "Residual",
      "Strongly_resolving" = "Strongly resolving"
    ))
  ) +
  
  scale_x_continuous(
    breaks = c(
      4,
      7,
      11,
      23
    )
  ) +
  
  labs(
    x = "Day",
    y = "Paclitaxel vs Vehicle log2 fold change",
    title =
      "Longitudinal transcriptional trajectories in DRG"
  ) +
  
  theme_classic()


trajectory_plot


# ---------------------------------------------------------------------
# 36. Save trajectory figure
# ---------------------------------------------------------------------

ggsave(
  filename =
    "figures/GSE185084_DRG_gene_trajectory_classes.pdf",
  plot = trajectory_plot,
  width = 9,
  height = 6
)


# =====================================================================
# PART I
# FINAL VALIDATION CHECKS
# =====================================================================


cat(
  "\nDay-specific significant genes:\n"
)

print(
  c(
    Day4 = sig_day4,
    Day7 = sig_day7,
    Day11 = sig_day11,
    Day23 = sig_day23
  )
)


cat(
  "\nSignificant treatment x time interaction genes:\n"
)

print(
  interaction_count
)


cat(
  "\nFinal Day-11 trajectory classes:\n"
)

print(
  table(
    trajectory$final_class
  )
)


cat(
  "\nResolving genes with sign-crossing overshoot:\n"
)

print(
  sum(
    trajectory$resolving_overshoot
  )
)


cat(
  "\nScript 06 completed successfully.\n"
)