# Define function to return top and bottom genes of a DESeqResults object by LFC
# Or return all genes if asked.

extract_topgenes <- function(results, dds, ntop = Inf, signif_only = TRUE) {
  # Format results and gene metadata
  results_format <- results %>%
    as_tibble(rownames = "gene_id") %>%
    dplyr::select(c(gene_id, log2FoldChange, padj))

  row_ranges_format <- rowRanges(dds) %>%
    as_tibble() %>%
    dplyr::select(c(gene_id, gene_name, description))

  # Merge into one table, sort by decreasing fold change
  results_merge <- right_join(
    row_ranges_format, results_format,
    by = join_by(gene_id == gene_id)
  ) %>%
    drop_na() %>%
    arrange(desc(log2FoldChange))

  # Remove non-significant genes
  if (signif_only == TRUE) {
    results_merge %<>% filter(padj < 0.05)
  }

  results_merge %<>% dplyr::rename(
    "ENSEMBL ID" = gene_id,
    "Symbol" = gene_name,
    "Gene name" = description,
    "LogFC" = log2FoldChange,
    "adj. P-val" = padj
  )

  # Return split table (or not) depending on number of top genes provided
  # If not provided: Return entire table
  # If smaller number of top genes, return table containing top and bottom genes

  if (ntop == Inf) {
    results_merge %>%
      return()
  } else {
    results_merge %$%
      rbind(
        filter(., `LogFC` > 0) %>% slice_head(n = ntop),
        filter(., `LogFC` < 0) %>% slice_tail(n = ntop)
      ) %>%
      return()
  }
}

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
