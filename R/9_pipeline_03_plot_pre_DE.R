here::i_am("R/9_pipeline_03_plot_pre_DE.R")

########################
# Run figures for pre-DE
########################

# Run scripts of pipeline
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
      source(here::here("R", x), echo = TRUE, verbose = FALSE)
    },
    .progress = TRUE
  )
