## Implement MA plots 2D volcano plots that support outliers
# Define function to clip logFC; padj to desired point
# And clip padj
# Add relevant metadata

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
    .default = 2
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

# Hardcode plot specs

cutoff_logfc <- 3
cutoff_padj <- 1e-15
alpha <- 0.05

# Load data

quant_deseq2_batchcor <- readRDS(file = here::here(
  "output",
  "data_expression",
  "post_DGE",
  "quant_deseq2_batchcor.RDS"
))

results_lfcshrink <- readRDS(file = here::here(
  "output",
  "data_expression",
  "post_DGE",
  "results_deseq2_lfcshrink.RDS"
))

label_genes <- rowRanges(quant_deseq2_batchcor)$gene_name

# Clip data

results_clipped <- map(
  results_lfcshrink,
  \(x) {
    clip_results(
      x,
      cutoff_logfc = cutoff_logfc,
      cutoff_padj = cutoff_padj,
      alpha = alpha
    )
  },
  .progress = TRUE
)

# Draw 2D volcano plots

plots_volcano_ma <- purrr::imap(
  results_clipped,
  \(x, idx) {
    plot_grid(
      EnhancedVolcano(
        x,
        lab = label_genes,
        x = "log2FoldChange",
        y = "padj",
        xlim = c(cutoff_logfc * -1, cutoff_logfc),
        ylim = c(0, -log10(cutoff_padj)),
        ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
        axisLabSize = 12,
        title = idx,
        titleLabSize = 12,
        subtitle = NULL,
        caption = NULL,
        pCutoff = 0.05,
        FCcutoff = 0.5,
        pointSize = x$volcano_size,
        labSize = 2,
        boxedLabels = FALSE,
        shapeCustom = x$volcano_shape,
        legendPosition = "none",
        drawConnectors = TRUE,
        widthConnectors = 0.1,
        colConnectors = "grey60",
        arrowheads = FALSE
      ),
      EnhancedVolcano(
        x,
        lab = label_genes,
        x = "log2FoldChange",
        y = "baseMean_new",
        selectLab = label_genes[which(x$padj <= alpha)],
        xlim = c(cutoff_logfc * -1, cutoff_logfc),
        ylim = c(0.5, 7),
        ylab = bquote(~ -Log[10] ~ "mean of normalised counts"),
        axisLabSize = 12,
        title = idx,
        titleLabSize = 12,
        subtitle = NULL,
        caption = NULL,
        pCutoff = NA,
        FCcutoff = NA,
        pointSize = x$volcano_size,
        labSize = 2,
        boxedLabels = FALSE,
        colCustom = x$colour,
        shapeCustom = x$volcano_shape,
        legendPosition = "none",
        drawConnectors = TRUE,
        widthConnectors = 0.1,
        colConnectors = "grey60",
        arrowheads = FALSE
      ),
      ncol = 2,
      nrow = 1
    )
  },
  .progress = TRUE
)
