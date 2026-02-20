here::i_am("R/0_define_results_extract.R")

########################
# Define function for extraction of top genes from a DESeqResults object
# sorted by LFC
# Can select top and bottom genes, or simply return all genes.
########################

# Import packages
library(conflicted)
library(DESeq2)
library(here)
library(tidyverse)

# Define function
extract_topgenes <- function(
  results,
  dds,
  ntop = Inf,
  signif_only = TRUE
) {
  # Format results and gene metadata
  results_format <- format_deseq_results(results)
  rowranges_format <- format_deseq_rowranges(dds)

  # Merge into one table, sort by decreasing fold change
  results_merge <- dplyr::right_join(
    x = rowranges_format,
    y = results_format,
    by = dplyr::join_by(gene_id == gene_id)
  ) %>%
    tidyr::drop_na() %>%
    dplyr::arrange(dplyr::desc(log2FoldChange))

  # Remove non-significant genes
  if (signif_only == TRUE) {
    results_merge <- results_merge %>%
      dplyr::filter(padj < 0.05)
  }

  # Rename columns to be human readable
  results_merge <- results_merge %>%
    dplyr::rename(
      "ENSEMBL ID" = gene_id,
      "Symbol" = gene_name,
      "Gene name" = description,
      "ENTREZ ID" = entrezid,
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
        dplyr::filter(., `LogFC` > 0) %>%
          dplyr::slice_head(n = ntop),
        dplyr::filter(., `LogFC` < 0) %>%
          dplyr::slice_tail(n = ntop)
      ) %>%
      return()
  }
}

# Define function to join 2 DESeqResults from the list of DESeqResults
extract_joined_results <- function(
  results_1,
  results_2,
  results_lrt,
  name_1,
  name_2,
  dds
) {
  # Format results and gene metadata
  rowranges_format <- format_deseq_rowranges(dds)

  results_1_format <- format_deseq_results(results_1)
  results_2_format <- format_deseq_results(results_2)

  results_lrt_format <- format_deseq_results(results_lrt) %>%
    dplyr::select(c(gene_id, padj)) %>%
    dplyr::rename(padj_lrt = padj)

  # Construct suffix for merged results
  suffix <- c(name_1, name_2) %>%
    purrr::map_chr(\(x) stringr::str_c("_", x))

  # Construct merged results table with inner join
  results_merge <- dplyr::inner_join(
    x = results_1_format,
    y = results_2_format,
    by = dplyr::join_by(gene_id == gene_id),
    suffix = suffix
  ) %>%
    dplyr::right_join(
      x = results_lrt_format,
      y = .,
      by = dplyr::join_by(gene_id == gene_id)
    ) %>%
    dplyr::right_join(
      x = rowranges_format,
      y = .,
      by = dplyr::join_by(gene_id == gene_id)
    ) %>%
    dplyr::rename(
      "ENSEMBL ID" = gene_id,
      "Symbol" = gene_name,
      "Gene name" = description,
      "ENTREZ ID" = entrezid,
    )

  # Return data
  return(results_merge)
}

## Create function to join 3 DESeqResults from the list of DESeqResults
extract_joined_results_trio <- function(
  results_1,
  results_2,
  results_3,
  results_lrt,
  name_1,
  name_2,
  name_3,
  dds
) {
  # Define colnames to be changed throughout function
  cols_to_rename <- c("log2FoldChange", "padj")

  # Format results and gene metadata
  rowranges_format <- format_deseq_rowranges(dds)

  results_1_format <- format_deseq_results(results_1) %>%
    dplyr::rename_with(
      ~ str_c(.x, name_1, sep = "_"),
      dplyr::all_of(cols_to_rename)
    )

  results_2_format <- format_deseq_results(results_2) %>%
    dplyr::rename_with(
      ~ str_c(.x, name_2, sep = "_"),
      dplyr::all_of(cols_to_rename)
    )

  results_3_format <- format_deseq_results(results_3) %>%
    dplyr::rename_with(
      ~ str_c(.x, name_3, sep = "_"),
      dplyr::all_of(cols_to_rename)
    )

  results_lrt_format <- format_deseq_results(results_lrt) %>%
    dplyr::select(c(gene_id, padj)) %>%
    dplyr::rename(padj_lrt = padj)

  # Construct merged results table with inner join
  results_merge <- dplyr::inner_join(
    x = results_1_format,
    y = results_2_format,
    by = dplyr::join_by(gene_id == gene_id)
  ) %>%
    dplyr::inner_join(
      x = .,
      y = results_3_format,
      by = join_by(gene_id == gene_id)
    ) %>%
    dplyr::right_join(
      x = results_lrt_format,
      y = .,
      by = join_by(gene_id == gene_id)
    ) %>%
    dplyr::right_join(
      x = rowranges_format,
      y = .,
      by = join_by(gene_id == gene_id)
    ) %>%
    dplyr::rename(
      "ENSEMBL ID" = gene_id,
      "Symbol" = gene_name,
      "Gene name" = description,
      "ENTREZ ID" = entrezid,
    )

  # Return data
  return(results_merge)
}
