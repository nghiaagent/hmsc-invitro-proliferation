### Don't source this file by itself; call in from another file after running env_prep.R

# Load dataset, run DESeq2 analysis
# Remember that experimental design is already embedded in the dataset.

quant_deseq2 <- readRDS(file = "./output/data_expression/pre_DGE/quant_refseq_deseq.RDS")

# Run DESeq2

quant_deseq2 %<>% DESeq()

x <- quant_deseq2

# Save data

saveRDS(quant_deseq2,
        "./output/data_expression/pre_DGE/quant_refseq_deseq_analysed.RDS")