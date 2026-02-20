here::i_am("R/9_pipeline_02_dge.R")

########################
# Run DGE
########################

# Run scripts of DESeq2 pipeline
c(
  "2_dge.R",
  "2_dge_nofilt.R"
) %>%
  walk(
    \(x) {
      message(paste0(
        "Sourcing ",
        here::here("scripts", x)
      ))
      source(here::here("R", x), echo = TRUE, verbose = FALSE)
    },
    .progress = TRUE
  )
