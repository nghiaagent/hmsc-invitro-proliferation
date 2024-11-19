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

design(quant_deseq2_batchcor) <- ~ condition_ID + cell_line + run_date

counts(quant_deseq2_batchcor) <- quant_deseq2_batchcor %$%
  sva::ComBat_seq(
    counts(.),
    batch = colData(.)$run_date,
    covar_mod = model.matrix(~ condition_ID + cell_line, data = colData(.))
  ) %>%
  `storage.mode<-`(., "integer")

# Export rlog

rlog_deseq2_batchcor <- rlog(quant_deseq2_batchcor)

# Run DESeq2
## Run once with full model
## Run again with reduced model and LRT
## to get ANOVA-like statistic for 3D volcano plot
## (i.e. chance of gene differentially expressed at all)

quant_deseq2_batchcor <- quant_deseq2_batchcor %>%
  DESeq()

quant_deseq2_lrt <- quant_deseq2_batchcor %>%
  DESeq(
    test = "LRT",
    reduced = ~ cell_line + run_date
  )

# Obtain results

results <- map(
  list_contrasts_deseq2,
  \(x) {
    results(
      quant_deseq2_batchcor,
      contrast = x,
      filterFun = ihw,
      alpha = 0.05
    )
  },
  .progress = TRUE
)

results_lfcshrink <- map2(
  .x = list_contrasts_deseq2,
  .y = results,
  \(x, y) {
    lfcShrink(
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

fwrite(
  counts(quant_deseq2_batchcor),
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
  rlog_deseq2_batchcor,
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
