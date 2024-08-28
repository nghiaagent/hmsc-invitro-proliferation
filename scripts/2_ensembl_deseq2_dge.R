### Don't source this file by itself; call in from another file after running env_prep.R

# Load dataset
# Remember that experimental design is already embedded in the dataset and model - only need to extract comparisons

quant_deseq2 <- readRDS("output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")

# Try batch correction

quant_deseq2_batchcor <- quant_deseq2

design(quant_deseq2_batchcor) <- ~ condition_ID + cell_line

counts(quant_deseq2_batchcor) <- quant_deseq2_batchcor %$%
  sva::ComBat_seq(
    counts(.),
    batch = colData(.)$run_date
    ## Uncomment the below to preserve covariates in batch correction.
    ## Currently this produces NAs due to unbalanced design.
    # ,
    # covar_mod = model.matrix( ~ condition_ID + cell_line, data = colData(.))
  ) %>%
  `storage.mode<-`(., "integer")

# Run DESeq2 on datasets

quant_deseq2_batchcor %<>% DESeq()
quant_deseq2 %<>% DESeq()

# Obtain results

results <- map(
  list_contrasts_deseq2[c(1,3,5,13,14,15)],
  \ (x) results(
    quant_deseq2,
    contrast = x,
    filterFun = ihw,
    alpha = 0.05
  )
)

results_batchcor <- map(
  list_contrasts_deseq2[c(1,3,5,13,14,15)],
  \ (x) results(
    quant_deseq2_batchcor,
    contrast = x,
    filterFun = ihw,
    alpha = 0.05
  )
)