here::i_am("R/9_pipeline_DGE.R")

########################
# Run DGE
########################

# Run scripts of DESeq2 pipeline
c(
  "3_dge_extract_genes_ORA.R",
  "3_plot_heatmap.R",
  "3_plot_PCA.R"
) %>%
  walk(
    \(x) {
      message(paste0(
        "Sourcing ",
        here::here("scripts", x)
      ))
      source(here::here("R", x), echo = TRUE, verbose = TRUE)
    },
    .progress = TRUE
  )
