# Load data

quant_deseq2 <- readRDS(
  file = here::here(
    "output", "data_expression", "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

quant_deseq2_lrt <- readRDS(
  file = here::here(
    "output", "data_expression", "post_DGE",
    "quant_deseq2_LRT.RDS"
  )
)

results_lfcshrink <- readRDS(
  file = here::here(
    "output", "data_expression", "post_DGE",
    "results_deseq2_lfcshrink.RDS"
  )
)

# Define list of results to be exported

list_sheets <- results_lfcshrink %$%
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
    "Treatment @ P+13, D5" = .[[6]],
    "D5 - D3, P+5 Ctl" = .[[7]],
    "D5 - D3, P+5 Hep" = .[[10]],
    "D5 - D3, P+7 Ctl" = .[[8]],
    "D5 - D3, P+7 Hep" = .[[11]],
    "D5 - D3, P+13 Ctl" = .[[9]],
    "D5 - D3, P+13 Hep" = .[[12]]
  )

# Extract logFC and p-vals for comparisons between passages at D3 and D5
# Compile to list, coerce to data.frame, export to excel
# 2 lists: Contrasts and LRT

## All contrasts (Wald test)

write.xlsx(
  x = map(
    list_sheets,
    \(x) {
      extract_topgenes(
        results = x,
        dds = quant_deseq2,
        ntop = Inf,
        signif_only = FALSE
      )
    },
    .progress = TRUE
  ),
  file = here::here(
    "output",
    "data_expression",
    "genes_all_contrasts.xlsx"
  ),
  asTable = TRUE
)

## Omnibus (LRT)

write.xlsx(
  x = results(quant_deseq2_lrt) %>%
    extract_topgenes(
      results = .,
      dds = quant_deseq2,
      ntop = Inf,
      signif_only = FALSE
    ) %>%
    dplyr::arrange(`adj. P-val`) %>%
    dplyr::select(!LogFC),
  file = here::here(
    "output",
    "data_expression",
    "genes_all_LRT.xlsx"
  ),
  asTable = TRUE
)

# Extract top and bottom 20 genes for desired comparisons
# Compile to list, export to excel
# 2 lists: all genes and DEGs only

map2(
  c(TRUE, FALSE),
  c("genes_top_DE.xlsx", "genes_top_LFC.xlsx"),
  \(signif_only, filename) {
    write.xlsx(
      x = map(
        list_sheets,
        \(x) {
          extract_topgenes(
            x,
            quant_deseq2,
            ntop = 20,
            signif_only = signif_only
          )
        },
        .progress = TRUE
      ),
      file = here::here("output", "data_expression", filename),
      asTable = TRUE
    )
  }
)
