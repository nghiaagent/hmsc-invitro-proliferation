# Extract genes based on their direction of regulation
# between P5 - P7 and P7 - P13
# (e.g. up - up, down - down, down - up, up - down)

results <- readRDS(file = here::here(
  "output",
  "data_expression",
  "post_DGE",
  "results_deseq2.RDS"
))

up_P5vsP7 <- subset(results[[13]], padj < 0.05) %>%
  subset(log2FoldChange > 0)

down_P5vsP7 <- subset(results[[13]], padj < 0.05) %>%
  subset(log2FoldChange < 0)

up_P7vsP13 <- subset(results[[14]], padj < 0.05) %>%
  subset(log2FoldChange > 0)

down_P7vsP13 <- subset(results[[14]], padj < 0.05) %>%
  subset(log2FoldChange < 0)

genes_up_up <-
  quant_DGE_voom$genes[which(tests[, 13] == 1 & tests[, 14] == 1), ]

genes_up_down <-
  quant_DGE_voom$genes[which(tests[, 13] == 1 & tests[, 14] == -1), ]

genes_down_up <-
  quant_DGE_voom$genes[which(tests[, 13] == -1 & tests[, 14] == 1), ]

genes_down_down <-
  quant_DGE_voom$genes[which(tests[, 13] == -1 & tests[, 14] == -1), ]

# Extract genes based on Hep response

# Export lists

fwrite(
  genes_up_up$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_up_ENSEMBL.txt",
  na = ""
)

fwrite(
  genes_up_up$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_up_GENENAME.txt",
  na = ""
)


fwrite(
  genes_up_down$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_down_ENSEMBL.txt",
  na = ""
)

fwrite(
  genes_up_down$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_down_GENENAME.txt",
  na = ""
)

fwrite(
  genes_up_down$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_down_ENSEMBL.txt",
  na = ""
)

fwrite(
  genes_up_down$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_down_GENENAME.txt",
  na = ""
)

fwrite(
  genes_down_down$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_down_down_ENSEMBL.txt",
  na = ""
)

fwrite(
  genes_down_down$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_down_down_GENENAME.txt",
  na = ""
)
