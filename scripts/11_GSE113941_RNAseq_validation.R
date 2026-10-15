list.files(
  "data_raw/GSE113941",
  full.names = FALSE
)
# ------------------------------------------------------------
# 2. Read processed TPM matrix
# ------------------------------------------------------------

tpm_file <- "data_raw/GSE113941/GSE113941_uniform_process_tpm_vals.txt.gz"

tpm113941 <- read.delim(
  gzfile(tpm_file),
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

dim(tpm113941)

names(tpm113941)

head(tpm113941)
dim(tpm113941)
names(tpm113941)
head(tpm113941)
# ------------------------------------------------------------
# 3. Select independent DRG input samples
# ------------------------------------------------------------

normal_cols <- c(
  "DRG_normal_input1_TPM",
  "DRG_normal_input2_TPM",
  "DRG_normal_input3_TPM",
  "DRG_normal_input4_TPM"
)

cipn_cols <- c(
  "DRG_CIPN_input5_TPM",
  "DRG_CIPN_input6_TPM",
  "MERGED_DRG_CIPN_input7_TPM",
  "DRG_CIPN_input8_TPM"
)

validation113941 <- tpm113941[
  ,
  c(
    "gene_short_name",
    normal_cols,
    cipn_cols
  )
]

dim(validation113941)

head(validation113941)
validation_matrix <- as.matrix(
  validation113941[, -1]
)

rownames(validation_matrix) <- validation113941$gene_short_name

storage.mode(validation_matrix) <- "numeric"

dim(validation_matrix)

summary(as.vector(validation_matrix))

round(
  cor(
    log2(validation_matrix + 1),
    method = "pearson"
  ),
  4
)
validation_matrix <- as.matrix(
  validation113941[, -1]
)

rownames(validation_matrix) <- validation113941$gene_short_name

storage.mode(validation_matrix) <- "numeric"

dim(validation_matrix)

summary(as.vector(validation_matrix))

round(
  cor(
    log2(validation_matrix + 1),
    method = "pearson"
  ),
  4
)
# ------------------------------------------------------------
# 4. PCA of independent GSE113941 validation samples
# ------------------------------------------------------------

log_validation113941 <- log2(validation_matrix + 1)

keep113941 <- rowSums(validation_matrix >= 1) >= 2

log_validation113941_filt <- log_validation113941[
  keep113941,
  ,
  drop = FALSE
]

dim(log_validation113941_filt)

gene_var113941 <- apply(
  log_validation113941_filt,
  1,
  var
)

top500_113941 <- names(
  sort(
    gene_var113941,
    decreasing = TRUE
  )
)[1:500]

pca113941 <- prcomp(
  t(log_validation113941_filt[top500_113941, ]),
  scale. = FALSE
)

pca113941_percent <- round(
  100 * pca113941$sdev^2 /
    sum(pca113941$sdev^2),
  1
)

group113941 <- factor(
  c(
    rep("Normal", 4),
    rep("Paclitaxel", 4)
  ),
  levels = c("Normal", "Paclitaxel")
)

pca113941_df <- data.frame(
  sample = colnames(log_validation113941_filt),
  group = group113941,
  PC1 = pca113941$x[, 1],
  PC2 = pca113941$x[, 2]
)

pca113941_percent[1:2]

pca113941_df
library(ggplot2)

p_pca113941 <- ggplot(
  pca113941_df,
  aes(
    x = PC1,
    y = PC2,
    shape = group
  )
) +
  geom_point(size = 4) +
  geom_text(
    aes(label = sample),
    vjust = -0.8,
    size = 3
  ) +
  labs(
    x = paste0("PC1 (", pca113941_percent[1], "%)"),
    y = paste0("PC2 (", pca113941_percent[2], "%)"),
    title = "GSE113941 DRG validation PCA",
    subtitle = "Normal vs paclitaxel, top 500 variable genes"
  ) +
  theme_bw()

p_pca113941
pca113941_percent[1:2]

pca113941_df

p_pca113941
# ------------------------------------------------------------
# 5. Sample distance QC
# ------------------------------------------------------------

sample_dist113941 <- dist(
  t(log_validation113941_filt)
)

as.matrix(sample_dist113941)
hc113941 <- hclust(sample_dist113941)

plot(
  hc113941,
  main = "GSE113941 sample clustering",
  xlab = "",
  sub = ""
)
ggsave(
  "figures/GSE113941_DRG_PCA_top500.pdf",
  p_pca113941,
  width = 8,
  height = 6
)
# ------------------------------------------------------------
# 6. Targeted validation of unresolved GSE185084 genes
# ------------------------------------------------------------

driver_genes_185084 <- c(
  "Pdia6",
  "Hspa5",
  "Dbp",
  "Tef",
  "Bhlhe40",
  "Hlf"
)

driver113941 <- validation113941[
  validation113941$gene_short_name %in% driver_genes_185084,
]

driver113941
driver113941$normal_mean <- rowMeans(
  driver113941[, normal_cols],
  na.rm = TRUE
)

driver113941$cipn_mean <- rowMeans(
  driver113941[, cipn_cols],
  na.rm = TRUE
)

driver113941$log2FC_validation <- log2(
  (driver113941$cipn_mean + 1) /
    (driver113941$normal_mean + 1)
)

driver113941[
  ,
  c(
    "gene_short_name",
    "normal_mean",
    "cipn_mean",
    "log2FC_validation"
  )
]
# ------------------------------------------------------------
# 7. Limma statistics for GSE113941
# ------------------------------------------------------------

library(limma)

group113941 <- factor(
  c(
    rep("Normal", 4),
    rep("Paclitaxel", 4)
  ),
  levels = c("Normal", "Paclitaxel")
)

design113941 <- model.matrix(
  ~ group113941
)

colnames(design113941)

fit113941 <- lmFit(
  log_validation113941_filt,
  design113941
)

fit113941 <- eBayes(
  fit113941,
  robust = TRUE
)

res113941 <- topTable(
  fit113941,
  coef = "group113941Paclitaxel",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

dim(res113941)

head(res113941)
driver_stats113941 <- res113941[
  rownames(res113941) %in% driver_genes_185084,
  c(
    "logFC",
    "AveExpr",
    "P.Value",
    "adj.P.Val"
  )
]

driver_stats113941
# ------------------------------------------------------------
# 8. Restore gene symbols to limma results
# ------------------------------------------------------------

res113941_unsorted <- topTable(
  fit113941,
  coef = "group113941Paclitaxel",
  number = Inf,
  adjust.method = "BH",
  sort.by = "none"
)

res113941_unsorted$gene_short_name <-
  rownames(log_validation113941_filt)

driver_stats113941 <- res113941_unsorted[
  res113941_unsorted$gene_short_name %in% driver_genes_185084,
  c(
    "gene_short_name",
    "logFC",
    "AveExpr",
    "P.Value",
    "adj.P.Val"
  )
]

driver_stats113941
# ------------------------------------------------------------
# 9. Module-level validation using CAMERA
# ------------------------------------------------------------

library(limma)

er_stress_genes <- c(
  "Pdia6",
  "Hspa5"
)

rhythmic_genes <- c(
  "Dbp",
  "Hlf",
  "Bhlhe40",
  "Tef"
)

gene_sets113941 <- list(
  ER_stress_UPR = which(
    rownames(log_validation113941_filt) %in% er_stress_genes
  ),
  
  Rhythmic_process = which(
    rownames(log_validation113941_filt) %in% rhythmic_genes
  )
)

sapply(
  gene_sets113941,
  length
)

camera113941 <- camera(
  log_validation113941_filt,
  index = gene_sets113941,
  design = design113941,
  contrast = 2
)

camera113941
sapply(gene_sets113941, length)

camera113941
# ------------------------------------------------------------
# 10. Prepare ranked gene list for Hallmark GSEA
# ------------------------------------------------------------

res113941_unsorted$gene_short_name <-
  rownames(log_validation113941_filt)

rank113941_df <- res113941_unsorted[
  !is.na(res113941_unsorted$gene_short_name) &
    res113941_unsorted$gene_short_name != "",
  c(
    "gene_short_name",
    "t"
  )
]

# Remove duplicate gene symbols by retaining the probe/gene entry
# with the largest absolute t-statistic
rank113941_df <- rank113941_df[
  order(abs(rank113941_df$t), decreasing = TRUE),
]

rank113941_df <- rank113941_df[
  !duplicated(rank113941_df$gene_short_name),
]

rank113941 <- rank113941_df$t
names(rank113941) <- rank113941_df$gene_short_name

rank113941 <- sort(
  rank113941,
  decreasing = TRUE
)

length(rank113941)

head(rank113941)
tail(rank113941)
# ------------------------------------------------------------
# 11. Hallmark gene sets
# ------------------------------------------------------------

library(fgsea)
library(msigdbr)

hallmark_mouse_113941 <- msigdbr(
  species = "Mus musculus",
  collection = "H"
)

hallmark_pathways_113941 <- split(
  hallmark_mouse_113941$gene_symbol,
  hallmark_mouse_113941$gs_name
)

length(hallmark_pathways_113941)
# ------------------------------------------------------------
# 12. Hallmark GSEA for GSE113941
# ------------------------------------------------------------

set.seed(123)

fgsea113941 <- fgsea(
  pathways = hallmark_pathways_113941,
  stats = rank113941,
  minSize = 10,
  maxSize = 500
)

fgsea113941 <- fgsea113941[
  order(fgsea113941$padj),
]

sum(
  fgsea113941$padj < 0.05,
  na.rm = TRUE
)

head(
  fgsea113941[
    ,
    c(
      "pathway",
      "NES",
      "pval",
      "padj",
      "size"
    )
  ],
  20
)
# ------------------------------------------------------------
# 13. Re-run fgsea with increased simple permutations
# ------------------------------------------------------------

set.seed(123)

fgsea113941 <- fgsea(
  pathways = hallmark_pathways_113941,
  stats = rank113941,
  minSize = 10,
  maxSize = 500,
  nPermSimple = 10000
)

fgsea113941 <- fgsea113941[
  order(fgsea113941$padj),
]

sum(
  fgsea113941$padj < 0.05,
  na.rm = TRUE
)

sum(is.na(fgsea113941$NES))
fgsea113941_save <- as.data.frame(fgsea113941)

fgsea113941_save$leadingEdge <- sapply(
  fgsea113941_save$leadingEdge,
  paste,
  collapse = ";"
)

write.csv(
  fgsea113941_save,
  "results/GSE113941_GSEA_Hallmark.csv",
  row.names = FALSE
)
# ------------------------------------------------------------
# 14. Compare GSE113941 with GSE185084 day 11
# ------------------------------------------------------------

gsea185084_day11 <- read.csv(
  "results/GSE185084_GSEA_Hallmark_day11.csv",
  stringsAsFactors = FALSE
)

comparison_hallmark <- merge(
  as.data.frame(fgsea113941)[
    ,
    c("pathway", "NES", "padj")
  ],
  gsea185084_day11[
    ,
    c("pathway", "NES", "padj")
  ],
  by = "pathway",
  suffixes = c(
    "_GSE113941",
    "_GSE185084_day11"
  )
)

dim(comparison_hallmark)
cor(
  comparison_hallmark$NES_GSE113941,
  comparison_hallmark$NES_GSE185084_day11,
  use = "complete.obs",
  method = "spearman"
)
cor(
  comparison_hallmark$NES_GSE113941,
  comparison_hallmark$NES_GSE185084_day11,
  use = "complete.obs",
  method = "pearson"
)
comparison_hallmark$sig_113941 <-
  comparison_hallmark$padj_GSE113941 < 0.05

comparison_hallmark$sig_185084 <-
  comparison_hallmark$padj_GSE185084_day11 < 0.05

comparison_hallmark$same_direction <-
  sign(comparison_hallmark$NES_GSE113941) ==
  sign(comparison_hallmark$NES_GSE185084_day11)

table(
  comparison_hallmark$sig_113941,
  comparison_hallmark$sig_185084
)

table(
  comparison_hallmark$same_direction
)
replicated_hallmark <- comparison_hallmark[
  comparison_hallmark$sig_113941 &
    comparison_hallmark$sig_185084,
]

replicated_hallmark[
  order(replicated_hallmark$padj_GSE113941),
]
table(
  replicated_hallmark$same_direction
)
# ------------------------------------------------------------
# 15. Hallmark NES concordance plot
# ------------------------------------------------------------

comparison_hallmark$replication_class <- "Other"

comparison_hallmark$replication_class[
  comparison_hallmark$sig_113941 &
    comparison_hallmark$sig_185084 &
    comparison_hallmark$same_direction
] <- "Significant_same_direction"

comparison_hallmark$replication_class[
  comparison_hallmark$sig_113941 &
    comparison_hallmark$sig_185084 &
    !comparison_hallmark$same_direction
] <- "Significant_opposite_direction"

comparison_hallmark$label <- ifelse(
  comparison_hallmark$sig_113941 &
    comparison_hallmark$sig_185084,
  gsub("HALLMARK_", "", comparison_hallmark$pathway),
  ""
)

library(ggplot2)
library(ggrepel)

p_validation_hallmark <- ggplot(
  comparison_hallmark,
  aes(
    x = NES_GSE185084_day11,
    y = NES_GSE113941
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    aes(shape = replication_class),
    size = 3
  ) +
  ggrepel::geom_text_repel(
    aes(label = label),
    size = 3,
    max.overlaps = Inf
  ) +
  labs(
    x = "GSE185084 Day 11 Hallmark NES",
    y = "GSE113941 Hallmark NES",
    title = "Independent pathway-level reproducibility",
    subtitle = "Paclitaxel-treated mouse DRG"
  ) +
  theme_bw()

p_validation_hallmark
ggsave(
  "figures/GSE113941_vs_GSE185084_Hallmark_validation.pdf",
  p_validation_hallmark,
  width = 8,
  height = 6
)
write.csv(
  comparison_hallmark,
  "results/GSE113941_vs_GSE185084_Hallmark_comparison.csv",
  row.names = FALSE
)
# ------------------------------------------------------------
# 16. Improved Hallmark validation figure
# ------------------------------------------------------------

comparison_hallmark$short_label <- ""

comparison_hallmark$short_label[
  comparison_hallmark$pathway == "HALLMARK_MYOGENESIS"
] <- "Myogenesis"

comparison_hallmark$short_label[
  comparison_hallmark$pathway ==
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
] <- "EMT"

comparison_hallmark$short_label[
  comparison_hallmark$pathway ==
    "HALLMARK_KRAS_SIGNALING_UP"
] <- "KRAS signaling up"

comparison_hallmark$short_label[
  comparison_hallmark$pathway ==
    "HALLMARK_ADIPOGENESIS"
] <- "Adipogenesis"

comparison_hallmark$short_label[
  comparison_hallmark$pathway ==
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
] <- "Oxidative phosphorylation"

p_validation_hallmark <- ggplot(
  comparison_hallmark,
  aes(
    x = NES_GSE185084_day11,
    y = NES_GSE113941
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    aes(shape = replication_class),
    size = 3
  ) +
  ggrepel::geom_text_repel(
    aes(label = short_label),
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.4,
    max.overlaps = Inf,
    min.segment.length = 0
  ) +
  scale_x_continuous(
    expand = expansion(mult = 0.18)
  ) +
  scale_y_continuous(
    expand = expansion(mult = 0.18)
  ) +
  labs(
    x = "GSE185084 Day 11 Hallmark NES",
    y = "GSE113941 Hallmark NES",
    title = "Independent pathway-level reproducibility",
    subtitle = "Paclitaxel-treated mouse DRG"
  ) +
  theme_bw()

p_validation_hallmark
ggsave(
  "figures/GSE113941_vs_GSE185084_Hallmark_validation.pdf",
  p_validation_hallmark,
  width = 8,
  height = 6
)
# ------------------------------------------------------------
# 17. Final publication-style validation figure
# ------------------------------------------------------------

label_data <- comparison_hallmark[
  comparison_hallmark$short_label != "",
]

label_data$label_x <- label_data$NES_GSE185084_day11
label_data$label_y <- label_data$NES_GSE113941

# Manual label positions
label_data$label_x[label_data$short_label == "Myogenesis"] <- 1.10
label_data$label_y[label_data$short_label == "Myogenesis"] <- 1.90

label_data$label_x[label_data$short_label == "EMT"] <- 1.55
label_data$label_y[label_data$short_label == "EMT"] <- 1.38

label_data$label_x[label_data$short_label == "KRAS signaling up"] <- -0.20
label_data$label_y[label_data$short_label == "KRAS signaling up"] <- 2.05

label_data$label_x[label_data$short_label == "Adipogenesis"] <- -2.30
label_data$label_y[label_data$short_label == "Adipogenesis"] <- 1.75

label_data$label_x[
  label_data$short_label == "Oxidative phosphorylation"
] <- -2.45

label_data$label_y[
  label_data$short_label == "Oxidative phosphorylation"
] <- 1.28


p_validation_hallmark_final <- ggplot(
  comparison_hallmark,
  aes(
    x = NES_GSE185084_day11,
    y = NES_GSE113941
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    aes(shape = replication_class),
    size = 3
  ) +
  geom_text(
    data = label_data,
    aes(
      x = label_x,
      y = label_y,
      label = short_label
    ),
    inherit.aes = FALSE,
    size = 3.5
  ) +
  scale_shape_manual(
    values = c(
      "Other" = 16,
      "Significant_opposite_direction" = 17,
      "Significant_same_direction" = 15
    ),
    labels = c(
      "Other",
      "Significant, opposite direction",
      "Significant, same direction"
    )
  ) +
  labs(
    x = "GSE185084 Day 11 Hallmark NES",
    y = "GSE113941 Hallmark NES",
    title = "Independent pathway-level reproducibility",
    subtitle = "Paclitaxel-treated mouse DRG",
    shape = NULL
  ) +
  theme_bw()

p_validation_hallmark_final
ggsave(
  "figures/GSE113941_vs_GSE185084_Hallmark_validation_FINAL.pdf",
  p_validation_hallmark_final,
  width = 8,
  height = 6
)
# ------------------------------------------------------------
# 18. Cleaner final Hallmark validation figure
# ------------------------------------------------------------

label_data <- comparison_hallmark[
  comparison_hallmark$short_label != "",
]

p_validation_hallmark_final2 <- ggplot(
  comparison_hallmark,
  aes(
    x = NES_GSE185084_day11,
    y = NES_GSE113941
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    aes(shape = replication_class),
    size = 3
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(label = short_label),
    size = 3.5,
    box.padding = 0.8,
    point.padding = 0.5,
    max.overlaps = Inf,
    min.segment.length = 0,
    seed = 123
  ) +
  scale_shape_manual(
    values = c(
      "Other" = 16,
      "Significant_opposite_direction" = 17,
      "Significant_same_direction" = 15
    ),
    labels = c(
      "Other",
      "Significant, opposite direction",
      "Significant, same direction"
    )
  ) +
  coord_cartesian(
    xlim = c(-3.4, 3.2),
    ylim = c(-1.4, 2.35),
    clip = "off"
  ) +
  labs(
    x = "GSE185084 Day 11 Hallmark NES",
    y = "GSE113941 Hallmark NES",
    title = "Independent pathway-level reproducibility",
    subtitle = "Paclitaxel-treated mouse DRG",
    shape = NULL
  ) +
  theme_bw() +
  theme(
    plot.margin = margin(10, 25, 10, 25)
  )

p_validation_hallmark_final2
ggsave(
  "figures/GSE113941_vs_GSE185084_Hallmark_validation_FINAL.pdf",
  p_validation_hallmark_final2,
  width = 8,
  height = 6
)
# ------------------------------------------------------------
# 19. Publication-style validation figure with fixed labels
# ------------------------------------------------------------

label_data <- comparison_hallmark[
  comparison_hallmark$short_label != "",
]

label_data$label_x <- c(
  2.75,   # Myogenesis
  2.75,   # EMT
  2.75,   # KRAS signaling up
  -3.15,   # Adipogenesis
  -3.15    # Oxidative phosphorylation
)

label_data$label_y <- c(
  1.95,
  1.65,
  2.25,
  1.75,
  1.35
)

p_validation_hallmark_final3 <- ggplot(
  comparison_hallmark,
  aes(
    x = NES_GSE185084_day11,
    y = NES_GSE113941
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    aes(shape = replication_class),
    size = 3
  ) +
  geom_segment(
    data = label_data,
    aes(
      x = NES_GSE185084_day11,
      y = NES_GSE113941,
      xend = label_x,
      yend = label_y
    ),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = label_data,
    aes(
      x = label_x,
      y = label_y,
      label = short_label
    ),
    inherit.aes = FALSE,
    hjust = ifelse(label_data$label_x > 0, 0, 1),
    size = 3.5
  ) +
  scale_shape_manual(
    values = c(
      "Other" = 16,
      "Significant_opposite_direction" = 17,
      "Significant_same_direction" = 15
    ),
    labels = c(
      "Other",
      "Significant, opposite direction",
      "Significant, same direction"
    )
  ) +
  coord_cartesian(
    xlim = c(-4, 4),
    ylim = c(-1.4, 2.5),
    clip = "off"
  ) +
  labs(
    x = "GSE185084 Day 11 Hallmark NES",
    y = "GSE113941 Hallmark NES",
    title = "Independent pathway-level reproducibility",
    subtitle = "Paclitaxel-treated mouse DRG",
    shape = NULL
  ) +
  theme_bw() +
  theme(
    plot.margin = margin(10, 40, 10, 40),
    legend.position = "right"
  )

p_validation_hallmark_final3
ggsave(
  "figures/GSE113941_vs_GSE185084_Hallmark_validation_FINAL3.pdf",
  p_validation_hallmark_final3,
  width = 9,
  height = 6
)
# ------------------------------------------------------------
# 20. Clean compact validation figure
# ------------------------------------------------------------

comparison_hallmark$plot_label <- ""

comparison_hallmark$plot_label[
  comparison_hallmark$pathway == "HALLMARK_MYOGENESIS"
] <- "Myogenesis"

comparison_hallmark$plot_label[
  comparison_hallmark$pathway ==
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
] <- "EMT"

comparison_hallmark$plot_label[
  comparison_hallmark$pathway ==
    "HALLMARK_KRAS_SIGNALING_UP"
] <- "KRAS-up"

comparison_hallmark$plot_label[
  comparison_hallmark$pathway ==
    "HALLMARK_ADIPOGENESIS"
] <- "Adipogenesis"

comparison_hallmark$plot_label[
  comparison_hallmark$pathway ==
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
] <- "OxPhos"

label_data <- comparison_hallmark[
  comparison_hallmark$plot_label != "",
]

p_validation_hallmark_clean <- ggplot(
  comparison_hallmark,
  aes(
    x = NES_GSE185084_day11,
    y = NES_GSE113941
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_point(
    aes(shape = replication_class),
    size = 3
  ) +
  ggrepel::geom_text_repel(
    data = label_data,
    aes(label = plot_label),
    size = 3.5,
    box.padding = 0.6,
    point.padding = 0.4,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  scale_shape_manual(
    values = c(
      "Other" = 16,
      "Significant_opposite_direction" = 17,
      "Significant_same_direction" = 15
    ),
    labels = c(
      "Other",
      "Significant, opposite direction",
      "Significant, same direction"
    )
  ) +
  labs(
    x = "GSE185084 Day 11 Hallmark NES",
    y = "GSE113941 Hallmark NES",
    title = "Independent pathway-level reproducibility",
    subtitle = "Paclitaxel-treated mouse DRG",
    shape = NULL
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(
    shape = guide_legend(nrow = 1)
  )

p_validation_hallmark_clean
ggsave(
  "figures/GSE113941_vs_GSE185084_Hallmark_validation_CLEAN.pdf",
  p_validation_hallmark_clean,
  width = 9,
  height = 6.5
)
# ------------------------------------------------------------
# 21. Save concise external validation summary
# ------------------------------------------------------------

validation_summary <- replicated_hallmark[
  ,
  c(
    "pathway",
    "NES_GSE113941",
    "padj_GSE113941",
    "NES_GSE185084_day11",
    "padj_GSE185084_day11",
    "same_direction"
  )
]

validation_summary$interpretation <- ifelse(
  validation_summary$same_direction,
  "Concordant significant reproducibility",
  "Significant in both, opposite direction"
)

validation_summary <- validation_summary[
  order(
    !validation_summary$same_direction,
    validation_summary$padj_GSE113941
  ),
]

validation_summary
write.csv(
  validation_summary,
  "results/GSE113941_external_validation_summary.csv",
  row.names = FALSE
)
# ------------------------------------------------------------
# 22. Save overall Hallmark validation metrics
# ------------------------------------------------------------

validation_metrics <- data.frame(
  metric = c(
    "Spearman_NES_correlation",
    "Pearson_NES_correlation",
    "Same_direction_pathways",
    "Opposite_direction_pathways",
    "Significant_GSE113941",
    "Significant_GSE185084_day11",
    "Significant_in_both",
    "Significant_same_direction",
    "Significant_opposite_direction"
  ),
  value = c(
    cor(
      comparison_hallmark$NES_GSE113941,
      comparison_hallmark$NES_GSE185084_day11,
      use = "complete.obs",
      method = "spearman"
    ),
    cor(
      comparison_hallmark$NES_GSE113941,
      comparison_hallmark$NES_GSE185084_day11,
      use = "complete.obs",
      method = "pearson"
    ),
    sum(comparison_hallmark$same_direction, na.rm = TRUE),
    sum(!comparison_hallmark$same_direction, na.rm = TRUE),
    sum(comparison_hallmark$sig_113941, na.rm = TRUE),
    sum(comparison_hallmark$sig_185084, na.rm = TRUE),
    nrow(replicated_hallmark),
    sum(replicated_hallmark$same_direction, na.rm = TRUE),
    sum(!replicated_hallmark$same_direction, na.rm = TRUE)
  )
)

validation_metrics
write.csv(
  validation_metrics,
  "results/GSE113941_external_validation_metrics.csv",
  row.names = FALSE
)

# End of Script 11
