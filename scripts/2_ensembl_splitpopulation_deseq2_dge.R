### Don't source this file by itself; call in from another file after running env_prep.R

# Load dataset, split into 2 populations, subset to D3 run DESeq2 analysis
# Remember that experimental design is already embedded in the dataset and model.

quant_big <- readRDS("./output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")

quant_cell_line <- c("hMSC-20176", "hMSC-21558") %>%
  set_names() %>%
  # Split dataset based on cell population
  map(.f = \(x) quant_big[, quant_big$cell_line == x & quant_big$timepoint_ID != "P13D5"]) %>%
  # Refactor batch and cell line factors, drop levels not present in the population
  map(.f = \(x) {
    colData(x)$cell_line %<>% droplevels()
    colData(x)$Passage %<>% droplevels()
    colData(x)$Day %<>% droplevels()
    colData(x)$Treatment %<>% droplevels()
    colData(x)$run_date %<>% droplevels()
    colData(x)$ID %<>% droplevels()
    colData(x)$timepoint_ID %<>% droplevels()
    colData(x)$condition_ID %<>% droplevels()
    return(x)
  }) %>%
  # Update design
  map(.f = \(x) {
    design(x) <- ~ condition_ID + run_date
    return(x)
    }) %>%
  # Run DESeq2
  map(.f = \(x) DESeq(x))

res <- results(
  quant_cell_line[[1]],
  contrast = c("condition_ID", "P7D3Untreated", "P5D3Untreated"),
  alpha = 0.05,
  filterFun = ihw
) %>%
  as.data.frame() %>%
  filter(padj < 0.05)

# Save data

saveRDS(
  quant_cell_line,
  "./output/data_expression/pre_DGE/quant_cDNA_deseq_splitpopulation_analysed.RDS"
)