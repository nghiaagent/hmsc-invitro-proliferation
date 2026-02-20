quant_deseq2_batchcor <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

conditions_interest <- c(
  "P5D5Untreated",
  "P13D5Untreated"
)

quant_small <- quant_deseq2_batchcor %$%
  .[, colData(.)$condition_ID %in% conditions_interest]

# quant_small <- quant_small %$%
#    .[, colData(.)$cell_line == "hMSC-21558"]

# Extract rownames of GOIs

genes_sel <- rowRanges(quant_small) %>%
  as.data.frame() %>%
  filter(gene_name %in% c("GADD45B", "JUNB", "NFKB1", "MAP3K7"))

gene_counts <- plotCounts(
  quant_small,
  gene = "ENSG00000109320",
  intgroup = c("Passage"),
  returnData = TRUE
)

ggplot(
  gene_counts,
  aes(
    x = Passage,
    y = count,
    color = Passage
  )
) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  scale_y_log10() +
  theme_minimal_hgrid() +
  theme_classic() +
  labs(
    x = "Passage",
    y = "Normalised counts"
  ) +
  scale_color_manual(values = c("#156082", "#196b24"))

ggsave(filename = "temp.png", width = 8, height = 6, scale = 0.6)
