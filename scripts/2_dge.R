# Load dataset
# Remember that experimental design is already embedded
# in the dataset and model - only need to extract comparisons

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

# Export rlog for camera; WGCNA, etc.

rlog_deseq2_batchcor <- rlog(quant_deseq2_batchcor)

# Save data

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
  rlog_deseq2_batchcor,
  file = here::here("output", "data_expression", "post_DGE", "rlog_deseq2.RDS")
)

saveRDS(
  results,
  file = here::here("output", "data_expression", "post_DGE", "results_deseq2.RDS")
)
