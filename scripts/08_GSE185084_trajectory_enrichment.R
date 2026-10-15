# ============================================================
# CIPN Computational Project
# Script 08: Functional enrichment of trajectory classes
# GSE185084 DRG
# ============================================================

# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(ggrepel)

# ------------------------------------------------------------
# 2. Load final longitudinal trajectory results
# ------------------------------------------------------------

trajectory <- read.csv(
  "results/GSE185084_longitudinal_gene_trajectories_FINAL.csv",
  stringsAsFactors = FALSE
)

dim(trajectory)
table(trajectory$final_class)

# ------------------------------------------------------------
# 3. Define trajectory gene sets
# ------------------------------------------------------------

genes_resolving <- trajectory$gene[
  trajectory$final_class == "Strongly_resolving"
]

genes_residual <- trajectory$gene[
  trajectory$final_class == "Residual"
]

genes_persistent <- trajectory$gene[
  trajectory$final_class == "Persistent"
]

genes_unresolved <- c(
  genes_persistent,
  genes_residual
)

length(genes_resolving)
length(genes_persistent)
length(genes_residual)
length(genes_unresolved)

head(genes_resolving)
genes_persistent
genes_residual

# ------------------------------------------------------------
# 4. Define enrichment background
# ------------------------------------------------------------

background_genes <- unique(trajectory$gene)

length(background_genes)

# ------------------------------------------------------------
# 5. Map SYMBOL to ENTREZID
# ------------------------------------------------------------

background_map <- bitr(
  background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

resolving_map <- bitr(
  genes_resolving,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

unresolved_map <- bitr(
  genes_unresolved,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

c(
  background_mapped =
    length(unique(background_map$ENTREZID)),
  resolving_mapped =
    length(unique(resolving_map$ENTREZID)),
  unresolved_mapped =
    length(unique(unresolved_map$ENTREZID))
)

setdiff(
  genes_resolving,
  resolving_map$SYMBOL
)

setdiff(
  genes_unresolved,
  unresolved_map$SYMBOL
)

# ------------------------------------------------------------
# 6. GO Biological Process enrichment
# ------------------------------------------------------------

ego_resolving <- enrichGO(
  gene = unique(resolving_map$ENTREZID),
  universe = unique(background_map$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

ego_unresolved <- enrichGO(
  gene = unique(unresolved_map$ENTREZID),
  universe = unique(background_map$ENTREZID),
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.20,
  readable = TRUE
)

dim(as.data.frame(ego_resolving))
dim(as.data.frame(ego_unresolved))

# ------------------------------------------------------------
# 7. Save GO enrichment results
# ------------------------------------------------------------

unresolved_GO <- as.data.frame(
  ego_unresolved
)

write.csv(
  unresolved_GO,
  "results/GSE185084_GO_BP_unresolved_genes.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_resolving),
  "results/GSE185084_GO_BP_resolving_genes.csv",
  row.names = FALSE
)

unresolved_GO[
  ,
  c(
    "Description",
    "geneID",
    "p.adjust"
  )
]

# ------------------------------------------------------------
# 8. Identify unresolved GO driver genes
# ------------------------------------------------------------

go_driver_genes <- unique(
  unlist(
    strsplit(
      unresolved_GO$geneID,
      "/"
    )
  )
)

go_driver_genes

trajectory[
  trajectory$gene %in% go_driver_genes,
  c(
    "gene",
    "final_class"
  )
]

# ------------------------------------------------------------
# 9. Summarize GO driver gene trajectories
# ------------------------------------------------------------

go_driver_summary <- trajectory[
  trajectory$gene %in% go_driver_genes,
  c(
    "gene",
    "log2FC_day11",
    "padj_day11",
    "log2FC_day23",
    "padj_day23",
    "relative_effect_day23",
    "final_class"
  )
]

go_driver_summary$functional_group <- ifelse(
  go_driver_summary$gene %in%
    c("Pdia6", "Hspa5"),
  "ER_proteostasis_associated",
  "Rhythmic_process"
)

go_driver_summary

write.csv(
  go_driver_summary,
  "results/GSE185084_unresolved_GO_driver_gene_trajectories.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Plot unresolved GO driver gene trajectories
# ------------------------------------------------------------

driver_plot_data <- data.frame(
  gene = rep(
    go_driver_summary$gene,
    each = 2
  ),
  day = rep(
    c("Day 11", "Day 23"),
    times = nrow(go_driver_summary)
  ),
  log2FC = as.vector(
    t(
      go_driver_summary[
        ,
        c(
          "log2FC_day11",
          "log2FC_day23"
        )
      ]
    )
  ),
  functional_group = rep(
    go_driver_summary$functional_group,
    each = 2
  )
)

driver_plot_data$day <- factor(
  driver_plot_data$day,
  levels = c(
    "Day 11",
    "Day 23"
  )
)

p_driver <- ggplot(
  driver_plot_data,
  aes(
    x = day,
    y = log2FC,
    group = gene
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  geom_line() +
  geom_point(
    size = 2.5
  ) +
  facet_wrap(
    ~ functional_group,
    scales = "free_y",
    labeller = as_labeller(c(
      "ER_proteostasis_associated" = "ER proteostasis-associated",
      "Rhythmic_process" = "Rhythmic-process associated"
    ))
  )  +
  ggrepel::geom_text_repel(
    aes(
      label = gene
    ),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.3,
    max.overlaps = Inf
  ) +
  labs(
    x = NULL,
    y = "Paclitaxel vs Vehicle log2 fold change",
    title = "Unresolved GO driver genes"
  ) +
  theme_bw()

p_driver

ggsave(
  "figures/GSE185084_unresolved_GO_driver_gene_trajectories.pdf",
  p_driver,
  width = 8,
  height = 5
)

# ------------------------------------------------------------
# 11. Final Script 08 summary
# ------------------------------------------------------------

cat(
  "\nScript 08 completed successfully.\n",
  "Strongly resolving genes:",
  length(genes_resolving),
  "\n",
  "Unresolved genes:",
  length(genes_unresolved),
  "\n",
  "Significant GO BP terms in resolving set:",
  nrow(as.data.frame(ego_resolving)),
  "\n",
  "Significant GO BP terms in unresolved set:",
  nrow(as.data.frame(ego_unresolved)),
  "\n"
)