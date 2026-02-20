here::i_am("R/2_dge_nofilt.R")

########################
# Perform DGE with minimal low-abundance filtering
# To try and keep GOIs that are lowly expressed for statistical analysis
########################

# Import packages
library(DESeq2)
library(edgeR)
library(here)
library(magrittr)
library(SummarizedExperiment)
library(tidyverse)

# Load dataset - no filter
quant_deseq2 <- readRDS(here::here(
  "output",
  "data_expression",
  "pre_DGE",
  "quant_cDNA_deseq_nofilter.RDS"
))

# Perform minimal filtering
keep <- quant_deseq2 %$%
  filterByExpr(
    y = counts(.),
    group = colData(.)$condition_ID,
    min.count = 10,
    min.total.count = 15
  )

quant_deseq2 <- quant_deseq2[keep, ]

# Perform batch correction
## Batch correction + factor
quant_deseq2_batchcor <- quant_deseq2
design(quant_deseq2_batchcor) <- ~ condition_ID + cell_line + run_date
counts(quant_deseq2_batchcor) <- quant_deseq2_batchcor %$%
  sva::ComBat_seq(
    counts(.),
    batch = colData(.)$run_date,
    covar_mod = model.matrix(~ condition_ID + cell_line, data = colData(.))
  ) %>%
  `storage.mode<-`(., "integer")

# Run DESeq2
quant_deseq2_batchcor <- quant_deseq2_batchcor %>%
  DESeq()

# Obtain results
results <- list_contrasts_deseq2 %>%
  map(
    \(contrast) {
      results(
        quant_deseq2_batchcor,
        contrast = contrast,
        alpha = 0.05
      )
    },
    .progress = TRUE
  )

# Save data
## DESeq dataset
saveRDS(
  quant_deseq2_batchcor,
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor_nofilter.RDS"
  )
)

## Results (non-shrink LFC)
saveRDS(
  results,
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2_nofilter.RDS"
  )
)
