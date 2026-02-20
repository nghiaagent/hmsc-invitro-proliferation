here::i_am("R/4_dge_plot_volcano3D.R")

########################
# Build 3D volcano plots for DGE
# Supply p-values to the 3D volcano plot
# First column: ANOVA p-values
# Remaining columns: Comparisons
# 13: A vs B (P5 vs P7)
# 15: A vs C (P5 vs P13)
# 14: B vs C (P7 vs P13)
########################

# Import packages
library(DESeq2)
library(here)
library(SummarizedExperiment)
library(tidyverse)

# Define axis specs
plotlist <- list(
  z_score = 1,
  logFC = 2
)

breaks <- seq(
  from = 0,
  to = 3.5,
  by = 0.5
)

# Load data
rlog_deseq2_batchcor <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "rlog_deseq2.RDS"
  )
)

results <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2.RDS"
  )
)

quant_deseq2_lrt <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_LRT.RDS"
  )
)

results_lrt <- results(quant_deseq2_lrt, filterFun = ihw, alpha = 0.05)

# Supply pvalues
polar_pvals <- cbind(
  results_lrt$pvalue,
  results[[13]]$pvalue,
  results[[15]]$pvalue,
  results[[14]]$pvalue
)

polar_padj <- cbind(
  results_lrt$padj,
  results[[13]]$padj,
  results[[15]]$padj,
  results[[14]]$padj
)

## Construct volcano3d object
outcome <- colData(rlog_deseq2_batchcor)$condition_ID %>%
  factor(
    levels = c(
      "P5D3Untreated",
      "P7D3Untreated",
      "P13D3Untreated"
    ),
    labels = c(
      "5",
      "7",
      "13"
    )
  )

polar_manual <- polar_coords(
  outcome = outcome,
  data = rlog_deseq2_batchcor %>%
    assay() %>%
    set_rownames(rowRanges(rlog_deseq2_batchcor)$gene_name) %>%
    t(),
  pvals = polar_pvals,
  padj = polar_padj,
  scheme = palette_volcano3d
)

rownames(polar_manual@pvals) <- rownames(rlog_deseq2_batchcor)
colnames(polar_manual@pvals) <- c("LRT", "P7vsP5", "P13vsP5", "P13vsP7")
rownames(polar_manual@padj) <- rownames(rlog_deseq2_batchcor)
colnames(polar_manual@padj) <- c("LRT", "P7vsP5", "P13vsP5", "P13vsP7")

## Plot
### Generate list of
### 3D volcano - Z-scores
### 3D volcano - logFC
### radial - Z-scores
### radial - logFC

volcano3d <- map(
  plotlist,
  \(x) {
    volcano3D(
      polar_manual,
      type = x,
      label_size = 24,
      z_axis_title_size = 20,
      radial_axis_title_size = 20
    )
  }
)

radial_plotly <- map(
  plotlist,
  \(x) {
    radial_plotly(
      polar_manual,
      type = x,
      r_axis_ticks = breaks
    )
  }
)

radial_ggplot <- map(
  plotlist,
  \(x) {
    radial_ggplot(
      polar_manual,
      type = x,
      r_axis_ticks = breaks
    )
  }
)

# Save data
saveRDS(
  polar_manual,
  file = here::here(
    "output",
    "plots_volcano3D",
    "polar_manual.RDS"
  )
)

saveRDS(
  list(
    volcano3d,
    radial_plotly,
    radial_ggplot
  ),
  file = here::here(
    "output",
    "plots_volcano3D",
    "volcano3d_passages_day3.RDS"
  )
)
