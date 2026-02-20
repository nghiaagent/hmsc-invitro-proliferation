here::i_am("R/0_define_results_clip.R")

########################
# Define function to clip logFC and padj in DESeqResults to desired point
# For plotting with EnhancedVolcano
# Add relevant metadata
########################

# Import packages
library(DESeq2)
library(here)
library(tidyverse)

# Define function
clip_results <- function(
  results,
  cutoff_logfc = 2.5,
  cutoff_padj = 1e-20,
  alpha = 0.05
) {
  # Clip logFC to the threshold
  results$log2FoldChange <- dplyr::case_when(
    results$log2FoldChange >= cutoff_logfc ~ cutoff_logfc,
    results$log2FoldChange <= (cutoff_logfc * -1) ~ (cutoff_logfc * -1),
    .default = results$log2FoldChange
  )

  # Clip padj to the threshold
  results$padj <- results$padj %>%
    dplyr::case_when(
      . <= cutoff_padj ~ cutoff_padj,
      .default = .
    )

  # Add custom shapes to dots to identify clipped genes
  ## Normal dots: Shape 19
  ## Clipped (positive): Shape -9658
  ## Clipped (negative): Shape -9668
  ## Add names for legend
  results$volcano_shape <- dplyr::case_when(
    results$log2FoldChange >= cutoff_logfc ~ -9658,
    results$log2FoldChange <= cutoff_logfc_neg ~ -9668,
    results$padj <= cutoff_padj ~ 17,
    .default = 19
  )

  names(results$volcano_shape) <- dplyr::case_when(
    results$volcano_shape == -9658 ~ str_c("logFC >", cutoff_logfc),
    results$volcano_shape == -9668 ~ str_c("logFC < -", cutoff_logfc),
    results$volcano_shape == 17 ~ str_c("padj <", cutoff_padj),
    results$volcano_shape == 19 ~ "Unclipped"
  )

  # Make clipped genes larger
  results$volcano_size <- dplyr::case_when(
    results$log2FoldChange >= cutoff_logfc ~ 3,
    results$log2FoldChange <= cutoff_logfc_neg ~ 3,
    results$padj <= cutoff_padj ~ 3,
    .default = 1
  )

  # Create new baseMean column for plotting with EnhancedVolcano
  results$baseMean_new <- 1 / (results$baseMean + 1)

  # Create new significance status level column
  results$colour <- dplyr::case_when(
    results$padj <= alpha ~ "red2",
    .default = "grey30"
  )

  names(results$colour) <- dplyr::case_when(
    results$colour == "red2" ~ str_c("padj <", alpha),
    results$colour == "grey30" ~ "ns"
  )

  # Return modified data
  return(results)
}
