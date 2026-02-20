here::i_am("R/9_pipeline_04_plot_post_DE.R")

########################
# Run plotting after DE
########################

# Run scripts of pipeline
c(
  "4_dge_plot_01_volcano3D.R",
  "4_dge_plot_02_volcano2d_and_ma.R",
  "4_dge_plot_03_quadrant.R",
  "4_dge_plot_04_GOI_boxplot.R",
  "4_dge_plot_05_GOI_boxplot_thesis_grid.R",
  "4_dge_table_topgenes_volcano3D.R",
  "4_dge_table_topgenes.R"
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
