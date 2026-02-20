here::i_am("R/3_dge_extract_genes_ORA.R")

########################
# Extract gene names from DESeq fit for ORA using g:Profiler
# All contrasts
########################

# Import packages
library(conflicted)
library(data.table)
library(DESeq2)
library(here)
library(tidyverse)

# Load data
quant_deseq2 <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

results_lfcshrink <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2_lfcshrink.RDS"
  )
)

# Format table object, get only gene ENSEMBL IDs, ENTREZ IDs, and symbol
# ENSEMBL IDs or ENTREZ IDs are to be entered into g:Profiler
# Symbol is to be entered into STRING
results_lfcshrink_format <- results_lfcshrink %>%
  purrr::map(
    \(x) {
      extract_topgenes(
        results = x,
        dds = quant_deseq2,
        ntop = Inf,
        signif_only = TRUE
      )
    },
    .progress = TRUE
  )

results_ensembl_ids <- results_lfcshrink_format %>%
  purrr::map(\(x) x$`ENSEMBL ID`)

results_entrez_ids <- results_lfcshrink_format %>%
  purrr::map(\(x) x$`ENTREZ ID`)

results_symbol <- results_lfcshrink_format %>%
  purrr::map(\(x) x$Symbol)

# Export data
## Export ENSEMBL IDs
purrr::imap(
  results_ensembl_ids,
  \(x, idx) {
    data.table::fwrite(
      list(x),
      file = here::here(
        "output",
        "list_genes",
        str_c(idx, "_ENSEMBL.txt")
      )
    )
  }
)

## Export gene names
purrr::imap(
  results_symbol,
  \(x, idx) {
    data.table::fwrite(
      list(x),
      file = here::here(
        "output",
        "list_genes",
        str_c(idx, "_symbol.txt")
      )
    )
  }
)

## Export ENTREZ ID
purrr::imap(
  results_entrez_ids,
  \(x, idx) {
    data.table::fwrite(
      list(x),
      file = here::here(
        "output",
        "list_genes",
        str_c(idx, "_ENTREZ.txt")
      )
    )
  }
)
