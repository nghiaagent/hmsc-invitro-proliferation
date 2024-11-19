## Define function to clip logFC and padj in DESeqResults to desired point
## For plotting with EnhancedVolcano
## Add relevant metadata
clip_results <- function(
    results,
    cutoff_logfc = 2.5,
    cutoff_padj = 1e-20,
    alpha = 0.05) {
    cutoff_logfc_neg <- cutoff_logfc * -1

    # Clip logFC to the threshold
    results$log2FoldChange %<>% case_when(
        . >= cutoff_logfc ~ cutoff_logfc,
        . <= cutoff_logfc_neg ~ cutoff_logfc_neg,
        .default = .
    )

    # Clip padj to the threshold
    results$padj %<>% case_when(
        . <= cutoff_padj ~ cutoff_padj,
        .default = .
    )

    # Add custom shapes to dots to identify clipped genes
    ## Normal dots: Shape 19
    ## Clipped (positive): Shape -9658
    ## Clipped (negative): Shape -9668
    ## Add names for legend
    results$volcano_shape <- case_when(
        results$log2FoldChange >= cutoff_logfc ~ -9658,
        results$log2FoldChange <= cutoff_logfc_neg ~ -9668,
        results$padj <= cutoff_padj ~ 17,
        .default = 19
    )

    names(results$volcano_shape) <- case_when(
        results$volcano_shape == -9658 ~ str_c("logFC >", cutoff_logfc),
        results$volcano_shape == -9668 ~ str_c("logFC < -", cutoff_logfc),
        results$volcano_shape == 17 ~ str_c("padj <", cutoff_padj),
        results$volcano_shape == 19 ~ "Unclipped"
    )

    # Make clipped genes larger
    results$volcano_size <- case_when(
        results$log2FoldChange >= cutoff_logfc ~ 3,
        results$log2FoldChange <= cutoff_logfc_neg ~ 3,
        results$padj <= cutoff_padj ~ 3,
        .default = 1
    )

    # Create new baseMean column for plotting with EnhancedVolcano
    results$baseMean_new <- 1 / (results$baseMean + 1)

    # Create new significance status level column
    results$colour <- case_when(
        results$padj <= alpha ~ "red2",
        .default = "grey30"
    )

    names(results$colour) <- case_when(
        results$colour == "red2" ~ str_c("padj <", alpha),
        results$colour == "grey30" ~ "ns"
    )

    # Return modified data
    return(results)
}
