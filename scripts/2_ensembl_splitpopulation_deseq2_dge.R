### Don't source this file by itself; call in from another file after running env_prep.R

# Load dataset, split into 2 populations, remove P13 D5 data to run DESeq2 analysis
# Remember that experimental design is already embedded in the dataset and model.
quant_cell_line <- list(
  "hMSC-20176" = "hMSC-20176",
  "hMSC-21558" = "hMSC-21558"
) %>%
  
  # Split dataset based on cell population
  map(.x = .,
      .f = \(x) {
    readRDS("./output/data_expression/pre_DGE/quant_cDNA_deseq.RDS") %>%
      .[, .$cell_line == x & .$timepoint_ID != "P13D5"]
    }) %>%
  
  # Refactor batch and cell line factors, drop levels not present in the population. Update design
  map(.f = \(x) {
    design(x) <- ~ condition_ID + run_date
    colData(x) %<>% droplevels()
    return(x)
  }) %>%
  
  # Update design
  map(.f = \(x) {
    return(x)
    }) %>%
  
  # Run DESeq
  map(.f = \(x) DESeq(x))

results_20176 <- map(
  list_contrasts_deseq2[c(1,3,5,13,14,15)],
  \ (x) results(
    quant_cell_line[["hMSC-20176"]],
    contrast = x,
    filterFun = ihw,
    alpha = 0.05
  ),
  .progress = TRUE
)

results_21558 <- map(
  list_contrasts_deseq2[c(1,3,5,13,14,15)],
  \ (x) results(
    quant_cell_line[["hMSC-21558"]],
    contrast = x,
    filterFun = ihw,
    alpha = 0.05
  ),
  .progress = TRUE
)

# Save data

saveRDS(
  quant_cell_line,
  "./output/data_expression/pre_DGE/quant_cDNA_deseq_splitpopulation_analysed.RDS"
)