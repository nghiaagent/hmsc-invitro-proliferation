here::i_am("R/2_dge.R")

########################
# Perform batch adjustment with ComBat-seq
# Then perform DESeq2 analysis
########################

# Import packages
library(ashr)
library(conflicted)
library(data.table)
library(DESeq2)
library(edgeR)
library(here)
library(IHW)
library(SummarizedExperiment)
library(tidyverse)

# Load dataset
# Remember that experimental design is already embedded
# in the dataset and model - only need to extract comparisons
quant_deseq2 <- readRDS(here::here(
  "output",
  "data_expression",
  "pre_DGE",
  "quant_cDNA_deseq.RDS"
))

# Perform batch correction
## Batch correction + factor
quant_deseq2_batchcor <- quant_deseq2
DESeq2::design(quant_deseq2_batchcor) <- ~ condition_ID + cell_line + run_date
DESeq2::counts(quant_deseq2_batchcor) <- quant_deseq2_batchcor %$%
  sva::ComBat_seq(
    DESeq2::counts(.),
    batch = SummarizedExperiment::colData(.)$run_date,
    covar_mod = model.matrix(~ condition_ID + cell_line, data = colData(.))
  ) %>%
  `storage.mode<-`(., "integer")

# Run DESeq2
## Run once with full model
## Run again with reduced model and LRT
## to get ANOVA-like statistic for 3D volcano plot
## (i.e. chance of gene differentially expressed at all)
## Run again with non-batch corrected data
quant_deseq2_batchcor <- quant_deseq2_batchcor %>%
  DESeq2::DESeq()

quant_deseq2_lrt <- quant_deseq2_batchcor %>%
  DESeq2::DESeq(
    test = "LRT",
    reduced = ~ cell_line + run_date
  )

quant_deseq2 <- quant_deseq2 %>%
  DESeq2::DESeq()

# Obtain results
## No LFC shrinking
results <- list_contrasts_deseq2 %>%
  purrr::map(
    \(x) {
      DESeq2::results(
        quant_deseq2_batchcor,
        contrast = x,
        filterFun = ihw,
        alpha = 0.05
      )
    },
    .progress = TRUE
  )

## With LFC shrinking
results_lfcshrink <- purrr::map2(
  .x = list_contrasts_deseq2,
  .y = results,
  \(x, y) {
    DESeq2::lfcShrink(
      quant_deseq2_batchcor,
      contrast = x,
      res = y,
      type = "ashr"
    )
  },
  .progress = TRUE
)

# Save data
## DESeq datasets
saveRDS(
  quant_deseq2,
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2.RDS"
  )
)

saveRDS(
  quant_deseq2_batchcor,
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

saveRDS(
  quant_deseq2_lrt,
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_LRT.RDS"
  )
)

### Export counts matrix as .csv
data.table::fwrite(
  DESeq2::counts(quant_deseq2_batchcor),
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor_counts.csv"
  ),
  row.names = TRUE,
  col.names = TRUE
)

## rlog
saveRDS(
  rlog(quant_deseq2_batchcor),
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "rlog_deseq2.RDS"
  )
)

## Results (non-shrink and shrunken LFC)
saveRDS(
  results,
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2.RDS"
  )
)

saveRDS(
  results_lfcshrink,
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2_lfcshrink.RDS"
  )
)
