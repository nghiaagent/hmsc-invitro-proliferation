here::i_am("R/0_define_results_format.R")

########################
# Define functions to format rowRanges and DESeqResults
# Used in other functions in the project
########################

# Import packages
library(DESeq2)
library(here)
library(tidyverse)

format_deseq_results <- function(results) {
  # Type check
  if (class(results) != "DESeqResults") {
    stop("Provided results must be DESeqResults class")
  }

  # Format data, extract only relevant columns
  results_format <- results %>%
    tibble::as_tibble(rownames = "gene_id") %>%
    dplyr::select(c(gene_id, log2FoldChange, padj))

  # Return data
  return(results_format)
}

format_deseq_rowranges <- function(dds) {
  # Type check
  if (class(dds) != "DESeqDataSet") {
    stop("Provided dataset must be DESeqDataSet class")
  }

  # Format data, extract only relevant columns
  rowranges_format <- rowRanges(dds) %>%
    tibble::as_tibble() %>%
    dplyr::select(c(gene_id, gene_name, entrezid, description, baseMean))

  # Return data
  return(rowranges_format)
}
