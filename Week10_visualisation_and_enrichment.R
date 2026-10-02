# Week 10: From a gene list to biology
# Input: DESeq2 results and normalised expression from Weeks 8 and 9
# Comparison: thin biofilm (45 days) vs early biofilm (10 days)
# Author: Kunjala Dissanayake
# Date: 18 September 2026

library(ggplot2)
library(pheatmap)
file.exists("Week9_results/deseq2_all_genes.csv")
file.exists("Week8_results/log_cpm.rds")
file.exists("Week8_results/sample_metadata.rds")
res_df <- read.csv("Week9_results/deseq2_all_genes.csv")
log_cpm <- readRDS("Week8_results/log_cpm.rds")
sample_metadata <- readRDS("Week8_results/sample_metadata.rds")
head(res_df)
dim(res_df)
str(res_df)
levels(sample_metadata$condition)

res_df$neg_log10_padj <- -log10(res_df$padj)

head(
  res_df[, c("gene_id", "log2FoldChange", "padj", "neg_log10_padj")]
)
res_df$status <- "Not significant"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up in thin"
res_df$status[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down in thin"

table(res_df$status)
p_volcano <- ggplot(
  res_df,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    colour = status
  )
) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_colour_manual(
    values = c(
      "Up in thin" = "#d73027",
      "Down in thin" = "#4575b4",
      "Not significant" = "grey70"
    )
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    colour = "grey40"
  ) +
  labs(
    title = "Differential expression: thin versus early biofilm",
    x = "log2 fold change (thin / early)",
    y = "-log10 adjusted p-value",
    colour = ""
  ) +
  theme_bw()

p_volcano
top10 <- res_df[order(res_df$padj), ][1:10, ]

p_volcano +
  geom_text(
    data = top10,
    aes(label = gene_id),
    colour = "black",
    size = 2.5,
    vjust = -0.8
  )
dir.create("Week10_results", showWarnings = FALSE)

ggsave(
  "Week10_results/volcano_plot.png",
  p_volcano,
  width = 7,
  height = 6,
  dpi = 300
)
res_df$status_strict <- "Not significant"

res_df$status_strict[
  res_df$padj < 0.01 & res_df$log2FoldChange > 1
] <- "Up in thin"

res_df$status_strict[
  res_df$padj < 0.01 & res_df$log2FoldChange < -1
] <- "Down in thin"

table(res_df$status_strict)
table(res_df$status)
sum(res_df$status_strict != res_df$status)
sig_ids <- res_df$gene_id[
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
]

sig_ids <- sig_ids[!is.na(sig_ids)]

length(sig_ids)

heat_matrix <- log_cpm[
  rownames(log_cpm) %in% sig_ids,
]

dim(heat_matrix)
annotation_col <- data.frame(
  Stage = sample_metadata$condition
)

rownames(annotation_col) <- sample_metadata$sample_id

annotation_col
identical(
  rownames(annotation_col),
  colnames(heat_matrix)
)
pheatmap(
  heat_matrix,
  scale = "row",
  show_rownames = FALSE,
  annotation_col = annotation_col,
  main = "Differentially expressed genes: thin vs early"
)
pheatmap(
  heat_matrix,
  scale = "none",
  show_rownames = FALSE,
  main = "Unscaled"
)
pheatmap(
  heat_matrix,
  scale = "row",
  show_rownames = FALSE,
  main = "Row scaled"
)
library(org.Sc.sgd.db)

columns(org.Sc.sgd.db)

res_df$gene_name <- mapIds(
  org.Sc.sgd.db,
  keys = res_df$gene_id,
  column = "COMMON",
  keytype = "ORF",
  multiVals = "first"
)

head(
  res_df[, c("gene_id", "gene_name", "log2FoldChange", "padj")],
  10
)
res_df$description <- mapIds(
  org.Sc.sgd.db,
  keys = res_df$gene_id,
  column = "DESCRIPTION",
  keytype = "ORF",
  multiVals = "first"
)

head(res_df$description, 3)
head(
  res_df[order(res_df$padj),
         c("gene_id", "gene_name", "log2FoldChange")],
  5
)
library(clusterProfiler)

gene_universe <- res_df$gene_id[!is.na(res_df$padj)]

length(gene_universe)
length(sig_ids)

ego <- enrichGO(
  gene = sig_ids,
  universe = gene_universe,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(
  as.data.frame(ego)[,
                     c("Description", "GeneRatio", "BgRatio", "p.adjust")
  ],
  10
)
dotplot(
  ego,
  showCategory = 15
) +
  labs(
    title = "Enriched biological processes: thin vs early"
  )
up_ids <- res_df$gene_id[
  res_df$padj < 0.05 &
    res_df$log2FoldChange > 1
]

up_ids <- up_ids[!is.na(up_ids)]

ego_up <- enrichGO(
  gene = up_ids,
  universe = gene_universe,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(as.data.frame(ego_up)$Description, 10)
down_ids <- res_df$gene_id[
  res_df$padj < 0.05 &
    res_df$log2FoldChange < -1
]

down_ids <- down_ids[!is.na(down_ids)]

ego_down <- enrichGO(
  gene = down_ids,
  universe = gene_universe,
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05
)

head(as.data.frame(ego_down)$Description, 10)
write.csv(
  res_df,
  "Week10_results/annotated_results.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego),
  "Week10_results/go_enrichment.csv",
  row.names = FALSE
)

saveRDS(
  heat_matrix,
  "Week10_results/heat_matrix.rds"
)

list.files("Week10_results")
top30 <- res_df[!is.na(res_df$padj), ]
top30 <- head(top30[order(top30$padj), ], 30)

top30_matrix <- log_cpm[
  top30$gene_id,
  ,
  drop = FALSE
]

top30_labels <- ifelse(
  is.na(top30$gene_name) | top30$gene_name == "",
  top30$gene_id,
  top30$gene_name
)

rownames(top30_matrix) <- make.unique(top30_labels)

pheatmap(
  top30_matrix,
  scale = "row",
  show_rownames = TRUE,
  annotation_col = annotation_col,
  main = "Top 30 differentially expressed genes"
)
png(
  "Week10_results/top30_heatmap.png",
  width = 8,
  height = 10,
  units = "in",
  res = 300
)

pheatmap(
  top30_matrix,
  scale = "row",
  show_rownames = TRUE,
  annotation_col = annotation_col,
  main = "Top 30 differentially expressed genes"
)

dev.off()
file.exists("Week10_results/top30_heatmap.png")
graphics.off()
dir.create("Week10_results", showWarnings = FALSE)

png(
  filename = "Week10_results/top30_heatmap.png",
  width = 2400,
  height = 3000,
  res = 300
)

pheatmap(
  top30_matrix,
  scale = "row",
  show_rownames = TRUE,
  annotation_col = annotation_col,
  main = "Top 30 differentially expressed genes"
)

dev.off()

file.exists("Week10_results/top30_heatmap.png")
up_results <- as.data.frame(ego_up)

up_top5 <- head(
  up_results[
    order(up_results$p.adjust),
    c("Description", "GeneRatio", "p.adjust")
  ],
  5
)

up_top5
write.csv(
  up_results,
  "Week10_results/go_enrichment_up.csv",
  row.names = FALSE
)
file.exists("Week10_results/volcano_plot.png")
file.exists("Week10_results/top30_heatmap.png")