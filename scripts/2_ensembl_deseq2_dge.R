### Don't source this file by itself; call in from another file after running env_prep.R

# Load dataset
# Remember that experimental design is already embedded in the dataset and model - only need to extract comparisons

quant_deseq2 <- readRDS("output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")

# Try batch correction

## Ignore batch correction

quant_deseq2_ignorebatch <- quant_deseq2

design(quant_deseq2_ignorebatch) <- ~ condition_ID + run_date

## Batch correction + factor

quant_deseq2_batchcor <- quant_deseq2

design(quant_deseq2_batchcor) <- ~ condition_ID + cell_line + run_date

counts(quant_deseq2_batchcor) <- quant_deseq2_batchcor %$%
  sva::ComBat_seq(
    counts(.),
    batch = colData(.)$run_date
    ## Uncomment the below to preserve covariates in batch correction.
    ## Currently this produces NAs due to unbalanced design.
    ,
    covar_mod = model.matrix( ~ condition_ID + cell_line, data = colData(.))
  ) %>%
  `storage.mode<-`(., "integer")

## Batch correction + no factor

quant_deseq2_batchcor_nofactor <- quant_deseq2

design(quant_deseq2_batchcor_nofactor) <- ~ condition_ID + cell_line

counts(quant_deseq2_batchcor_nofactor) <- quant_deseq2_batchcor_nofactor %$%
  sva::ComBat_seq(
    counts(.),
    batch = colData(.)$run_date
    ## Uncomment the below to preserve covariates in batch correction.
    ## Currently this produces NAs due to unbalanced design.
    ,
    covar_mod = model.matrix( ~ condition_ID + cell_line, data = colData(.))
  ) %>%
  `storage.mode<-`(., "integer")

## Batch correction + no covariate

quant_deseq2_batchcor_nocovar <- quant_deseq2

design(quant_deseq2_batchcor_nocovar) <- ~ condition_ID + cell_line + run_date

counts(quant_deseq2_batchcor_nocovar) <- quant_deseq2_batchcor_nocovar %$%
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

quant_deseq2 %<>% DESeq()
quant_deseq2_batchcor %<>% DESeq()
quant_deseq2_batchcor_nocovar %<>% DESeq()
quant_deseq2_batchcor_nofactor %<>% DESeq()
quant_deseq2_ignorebatch %<>% DESeq()

# Obtain results

## Compare between methods

results <- map(
  list(quant_deseq2,
       quant_deseq2_batchcor,
       quant_deseq2_batchcor_nocovar,
       quant_deseq2_batchcor_nofactor,
       quant_deseq2_ignorebatch
),
  \ (x) results(
    x,
    contrast = list_contrasts_deseq2[[15]],
    filterFun = ihw,
    alpha = 0.05
  ),
  .progress = TRUE
)

# 
# results <- map(
#   list_contrasts_deseq2[c(13,14,15)],
#   \ (x) results(
#     quant_deseq2,
#     contrast = x,
#     filterFun = ihw,
#     alpha = 0.05
#   ),
#   .progress = TRUE
# )
# 
# results_batchcor <- map(
#   list_contrasts_deseq2[c(13,14,15)],
#   \ (x) results(
#     quant_deseq2_batchcor,
#     contrast = x,
#     filterFun = ihw,
#     alpha = 0.05
#   ),
#   .progress = TRUE
# )