here::i_am("R/8_evaluate_batchcor.R")

########################
# Evaluate batch correction methods
########################

# Import packages
library(DESeq2)
library(here)
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
design(quant_deseq2_batchcor) <- ~ condition_ID + cell_line + run_date
counts(quant_deseq2_batchcor) <- quant_deseq2_batchcor %$%
  ComBat_seq(
    counts(.),
    batch = colData(.)$run_date,
    covar_mod = model.matrix(~ condition_ID + cell_line, data = colData(.))
  ) %>%
  `storage.mode<-`(., "integer")
# Create 3 datasets to evaluate batch effect correction method
## No batch correction, include batch term
## Batch correction and include batch term
## Batch correction no include batch term

## No batch correction, include batch term
quant_deseq2_nobatchcor_includeterm <- quant_deseq2 %>%
  DESeq()

## Batch correction, include batch term
quant_deseq2_batchcor_includeterm <- quant_deseq2_batchcor %>%
  DESeq()
## Batch correction, no batch term
quant_deseq2_batchcor_noterm <- quant_deseq2_batchcor
design(quant_deseq2_batchcor_noterm) <- ~ condition_ID + cell_line
quant_deseq2_batchcor_noterm <- quant_deseq2_batchcor_noterm %>%
  DESeq()

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
  imap(\(contrast, name_contrast) {
    imap(list_quant, \(quant, name_quant) {
      quant %>%
        results(
          contrast = contrast,
          filterFun = ihw,
          alpha = 0.05
        ) %>%
        as.data.frame() %>%
        rownames_to_column(var = "ensembl_id") %>%
        mutate(method = name_quant) %>%
        mutate(contrast = name_contrast)
    }) %>%
      rbindlist()
  }) %>%
  rbindlist() %>%
  # Convert method and contrasts to named factors for easier plotting
  mutate(
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
  ggplot(aes(x = x)) +
  geom_density(
    aes(
      x = pvalue,
      y = ..density..
    ),
    fill = "grey60",
    adjust = 0.5
  ) +
  geom_density(
    aes(
      x = padj,
      y = -..density..
    ),
    fill = "grey20",
    adjust = 0.5
  ) +
  geom_hline(yintercept = 0) +
  annotate(
    "label",
    x = 0.7,
    y = 5,
    label = "Raw p value",
    color = "grey40"
  ) +
  annotate(
    "label",
    x = 0.3,
    y = -10,
    label = "Adj. p value",
    color = "black"
  ) +
  labs(x = "P value") +
  theme_bw() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
  ) +
  facet_grid(contrast ~ method)

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
