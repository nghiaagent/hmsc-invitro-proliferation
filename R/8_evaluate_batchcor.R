here::i_am("R/8_evaluate_batchcor.R")

########################
# Evaluate batch correction methods
########################

# Import packages
library(conflicted)
library(DESeq2)
library(here)
library(IHW)
library(magrittr)
library(SummarizedExperiment)
library(sva)
library(tidyverse)

# Load data
quant_deseq2 <- readRDS(here::here(
  "output",
  "data_expression",
  "pre_DGE",
  "quant_cDNA_deseq.RDS"
))

# Apply batch correction to data
quant_deseq2_batchcor <- quant_deseq2
DESeq2::design(quant_deseq2_batchcor) <- ~ condition_ID + cell_line + run_date
DESeq2::counts(quant_deseq2_batchcor) <- quant_deseq2_batchcor %$%
  sva::ComBat_seq(
    DESeq2::counts(.),
    batch = SummarizedExperiment::colData(.)$run_date,
    covar_mod = model.matrix(~ condition_ID + cell_line, data = colData(.))
  ) %>%
  `storage.mode<-`(., "integer")
# Create 3 datasets to evaluate batch effect correction method
## No batch correction, include batch term
## Batch correction and include batch term
## Batch correction no include batch term

## No batch correction, include batch term
quant_deseq2_nobatchcor_includeterm <- quant_deseq2 %>%
  DESeq2::DESeq()

## Batch correction, include batch term
quant_deseq2_batchcor_includeterm <- quant_deseq2_batchcor %>%
  DESeq2::DESeq()
## Batch correction, no batch term
quant_deseq2_batchcor_noterm <- quant_deseq2_batchcor
DESeq2::design(quant_deseq2_batchcor_noterm) <- ~ condition_ID + cell_line
quant_deseq2_batchcor_noterm <- quant_deseq2_batchcor_noterm %>%
  DESeq2::DESeq()

# Create list of datasets for processing
list_quant <- list(
  nobatchcor_includeterm = quant_deseq2_nobatchcor_includeterm,
  batchcor_includeterm = quant_deseq2_batchcor_includeterm,
  batchcor_noterm = quant_deseq2_batchcor_noterm
)

# Get results of key contrasts
results_pilot <- list_contrasts_deseq2[
  c(
    "P7vsP5_UT_D3",
    "P13vsP7_UT_D3",
    "P13vsP5_UT_D3"
  )
] %>%
  purrr::imap(\(contrast, name_contrast) {
    purrr::imap(list_quant, \(quant, name_quant) {
      quant %>%
        DESeq2::results(
          contrast = contrast,
          filterFun = ihw,
          alpha = 0.05
        ) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "ensembl_id") %>%
        dplyr::mutate(method = name_quant) %>%
        dplyr::mutate(contrast = name_contrast)
    }) %>%
      data.table::rbindlist()
  }) %>%
  data.table::rbindlist() %>%
  # Convert method and contrasts to named factors for easier plotting
  dplyr::mutate(
    method = method %>%
      factor(
        levels = c(
          "nobatchcor_includeterm",
          "batchcor_noterm",
          "batchcor_includeterm"
        ),
        labels = c(
          "Include batch in model only",
          "Correct batch effect only",
          "Apply both"
        )
      ),
    contrast = contrast %>%
      factor(
        levels = c(
          "P7vsP5_UT_D3",
          "P13vsP7_UT_D3",
          "P13vsP5_UT_D3"
        ),
        labels = c(
          "P+7 vs. P+5",
          "P+13 vs. P+7",
          "P+13 vs. P+5"
        )
      )
  )

# Create density plot of raw p values
plot_pvals <- results_pilot %>%
  ggplot2::ggplot(ggplot2::aes(x = x)) +
  ggplot2::geom_density(
    ggplot2::aes(
      x = pvalue,
      y = ..density..
    ),
    fill = "grey60",
    adjust = 0.5
  ) +
  ggplot2::geom_density(
    ggplot2::aes(
      x = padj,
      y = -..density..
    ),
    fill = "grey20",
    adjust = 0.5
  ) +
  ggplot2::geom_hline(yintercept = 0) +
  ggpp::annotate(
    "label",
    x = 0.7,
    y = 5,
    label = "Raw p value",
    color = "grey40"
  ) +
  ggpp::annotate(
    "label",
    x = 0.3,
    y = -10,
    label = "Adj. p value",
    color = "black"
  ) +
  ggplot2::labs(x = "P value") +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.text.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
  ) +
  ggplot2::facet_grid(contrast ~ method)

# Save plots
## P value plot
ggsave(
  here::here(
    "output",
    "plots_QC",
    "Batch effect method comparison.png"
  ),
  plot_pvals,
  width = 8,
  height = 10,
  scale = 0.85
)
