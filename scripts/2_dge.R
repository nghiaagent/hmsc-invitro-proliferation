### Don't source this file by itself; call in from another file after running env_prep.R

# Load dataset
# Remember that experimental design is already embedded in the dataset and model - only need to extract comparisons

quant_deseq2 <- readRDS("output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")

# Perform batch correction
## Batch correction + factor

quant_deseq2_batchcor <- quant_deseq2

design(quant_deseq2_batchcor) <- ~ condition_ID + cell_line + run_date

counts(quant_deseq2_batchcor) <- quant_deseq2_batchcor %$%
  sva::ComBat_seq(
    counts(.),
    batch = colData(.)$run_date
    ,
    covar_mod = model.matrix( ~ condition_ID + cell_line, data = colData(.))
  ) %>%
  `storage.mode<-`(., "integer")

# Run DESeq2

quant_deseq2_batchcor %<>% DESeq()

# Obtain results

results <- map(
  list_contrasts_deseq2,
  \ (x) results(
    quant_deseq2_batchcor,
    contrast = x,
    filterFun = ihw,
    alpha = 0.05
  ),
  .progress = TRUE
)

results_batchcor <- map(
  list_contrasts_deseq2[c(1,3,5,13,14,15)],
  \ (x) results(
    quant_deseq2_batchcor,
    contrast = x,
    filterFun = ihw,
    alpha = 0.05
  ),
  .progress = TRUE
)