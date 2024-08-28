### Don't source this file by itself; call in from another file after running env_prep.R

# Load dataset, run DESeq2 analysis
# Remember that experimental design is already embedded in the dataset and model - only need to extract comparisons

quant_deseq2 <- readRDS("./output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")
quant_deseq2 %<>% DESeq()

quant_deseq2_batchcor <- readRDS("./output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")
design(quant_deseq2_batchcor) <-  ~ condition_ID + run_date + cell_line
counts(quant_deseq2_batchcor) <- quant_deseq2_batchcor %$%
  sva::ComBat_seq(counts(.), batch = colData(.)$run_date) %>%
  `storage.mode<-`(., "integer")
quant_deseq2_batchcor %<>% DESeq()


# Run DESeq2

results(quant_deseq2,
        contrast = c("condition_ID", "P5D3Untreated", "P5D3Treated"),
        filterFun = ihw,
        alpha = 0.05) %>% summary()

results(quant_deseq2,
        contrast = c("condition_ID", "P5D3Untreated", "P7D3Untreated"),
        filterFun = ihw,
        alpha = 0.05) %>% summary()

results(quant_deseq2,
        contrast = c("condition_ID", "P5D3Untreated", "P13D3Untreated"),
        filterFun = ihw,
        alpha = 0.05) %>% summary()

results(quant_deseq2_batchcor,
        contrast = c("condition_ID", "P5D3Untreated", "P5D3Treated"),
        filterFun = ihw,
        alpha = 0.05) %>% summary()

results(quant_deseq2_batchcor,
        contrast = c("condition_ID", "P5D3Untreated", "P7D3Untreated"),
        filterFun = ihw,
        alpha = 0.05) %>% summary()

results(quant_deseq2_batchcor,
        contrast = c("condition_ID", "P5D3Untreated", "P13D3Untreated"),
        filterFun = ihw,
        alpha = 0.05) %>% summary()


# Save data

saveRDS(quant_deseq2,
        "./output/data_expression/pre_DGE/quant_cDNA_deseq_analysed.RDS")