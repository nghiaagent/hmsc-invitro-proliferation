here::i_am("R/4_dge_plot_02_volcano2d_and_ma.R")

########################
# 2D volcano and MA plots
########################

# Import packages
library(conflicted)
library(cowplot)
library(DESeq2)
library(EnhancedVolcano)
library(here)
library(tidyverse)

# Hardcode plot specs
cutoff_logfc <- 3
cutoff_padj <- 1e-15
alpha <- 0.05

# Load data
quant_deseq2_batchcor <- readRDS(
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

label_genes <- SummarizedExperiment::rowRanges(quant_deseq2_batchcor)$gene_name

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
    cowplot::plot_grid(
      EnhancedVolcano::EnhancedVolcano(
        x,
        lab = label_genes,
        x = "log2FoldChange",
        y = "padj",
        xlim = c(cutoff_logfc * -1, cutoff_logfc),
        ylim = c(0, -log10(cutoff_padj)),
        ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
        axisLabSize = 8,
        title = idx,
        titleLabSize = 8,
        subtitle = NULL,
        caption = NULL,
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = x$volcano_size,
        labSize = 2,
        boxedLabels = FALSE,
        shapeCustom = x$volcano_shape,
        legendPosition = "none",
        drawConnectors = TRUE,
        widthConnectors = 0.1,
        colConnectors = "grey60",
        arrowheads = FALSE,
        gridlines.major = FALSE,
        gridlines.minor = FALSE
      ),
      EnhancedVolcano::EnhancedVolcano(
        x,
        lab = label_genes,
        x = "log2FoldChange",
        y = "baseMean_new",
        selectLab = label_genes[which(x$padj <= alpha)],
        xlim = c(cutoff_logfc * -1, cutoff_logfc),
        ylim = c(0.5, 7),
        ylab = bquote(~ Log[10] ~ "mean of normalised counts"),
        axisLabSize = 8,
        title = idx,
        titleLabSize = 8,
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
        arrowheads = FALSE,
        gridlines.major = FALSE,
        gridlines.minor = FALSE
      ),
      ncol = 2,
      nrow = 1,
      rel_widths = c(2, 1)
    )
  },
  .progress = TRUE
)

# Save data
saveRDS(
  plots_volcano_ma,
  here::here(
    "output",
    "plots_volcano2d",
    "plots_volcano2d_ma.RDS"
  )
)
