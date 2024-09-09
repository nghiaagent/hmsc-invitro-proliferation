### Don't source this file by itself; call in from another file after running env_prep.R

# Load dataset, run DESeq2 analysis
# Remember that experimental design is already embedded in the dataset and model - only need to extract comparisons

quant_deseq2 <- readRDS("./output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")

# Try subsetting dataset to just P5 D3

quant_subset <- quant_deseq2[, quant_deseq2$Passage == "P5" & quant_deseq2$Day == "D3" & quant_deseq2$cell_line == "hMSC-20176"]
design(quant_subset) <- ~ Treatment
quant_subset %<>% DESeq()

# Save data

saveRDS(quant_subset,
        "./output/data_expression/pre_DGE/quant_cDNA_deseq_subset_analysed.RDS")