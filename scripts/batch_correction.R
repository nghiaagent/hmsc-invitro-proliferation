# Apply ComBat-seq to control for batch effect
source("./scripts/dge.R")


quant_DGE_batchcor_withcovariates <-
  sva::ComBat_seq(quant_DGE_clean,
                  batch = quant_DGE_clean$samples$run_date,
                  covar_mod = design)

quant_DGE_batchcor_nocovariates <- sva::ComBat_seq(quant_DGE_clean$counts,
                                                   batch = table_design$run_date,
                                                   group = NULL,
                                                   full_mod = FALSE)
quant_DGE_clean_batchcor <- quant_DGE_clean
quant_DGE_clean_batchcor$counts <- quant_DGE_batchcor_nocovariates
