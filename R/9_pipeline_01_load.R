here::i_am("R/9_pipeline_01_load.R")

########################
# Run RNA-seq loading
########################

# Run scripts of pipeline
c(
  "1_pre_01_quant_import.R",
  "1_pre_02_QC_plots.R"
) %>%
  walk(
    \(x) {
      message(paste0("Sourcing ", here::here("scripts", x)))
      source(here::here("R", x), echo = TRUE, verbose = FALSE)
    },
    .progress = TRUE
  )
