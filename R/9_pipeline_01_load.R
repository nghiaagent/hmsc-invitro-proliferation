here::i_am("R/9_pipeline_01_load.R")

########################
# Run DGE
########################

# Run scripts of DESeq2 pipeline
c(
  "1_pre_01_quant_import.R",
  "1_pre_02_QC_plots.R"
) %>%
  walk(
    \(x) {
      message(paste0("Sourcing ", here::here("scripts", x)))
      source(here::here("R", x), echo = TRUE, verbose = TRUE)
    },
    .progress = TRUE
  )
