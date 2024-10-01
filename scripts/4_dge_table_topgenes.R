# Load data

results <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2.RDS"
  )
)

quant_deseq2 <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

# Extract logFC and p-vals for comparisons between passages at D3 and D5
# Review env prep file for info of which coefficient represent which contrasts
# Compile to list, coerce to data.frame, export to excel

gene_info <- rowRanges(quant_deseq2) %>%
  as.data.frame() %>%
  .[, c(6, 7, 8, 10)]

list_sheets <- results %$%
  list(
    "P+7 - P+5, D3" = .[[13]],
    "P+13 - P+5, D3" = .[[15]],
    "P+13 - P+7, D3" = .[[14]],
    "P+7 - P+5, D5" = .[[19]],
    "P+13 - P+5, D5" = .[[21]],
    "P+13 - P+7, D5" = .[[20]],
    "Treatment @ P+5, D3" = .[[1]],
    "Treatment @ P+7, D3" = .[[3]],
    "Treatment @ P+13, D3" = .[[5]],
    "Treatment @ P+5, D5" = .[[2]],
    "Treatment @ P+7, D5" = .[[4]],
    "Treatment @ P+13, D5" = .[[6]]
  ) %>%
  map(.f = \(x) {
    as.data.frame(x) %>%
      cbind(gene_info, .)
  })

write.xlsx(
  x = list_sheets,
  file = here::here(
    "output",
    "top_DEGs.xlsx"
  ),
  asTable = TRUE
)
