here::i_am("R/9_pipeline_07_GSVA.R")

########################
# Run GSVA analysis
########################

# Run scripts of pipeline
c(
  "7_post_GSVA_1.R",
  "7_post_GSVA_2_plot_volcano3D.R",
  "7_post_GSVA_3_resultstable.R",
  "7_post_GSVA_4_resultstable_volcano3d.R",
  "7_post_GSVA_5_POI_boxplot.R",
  "7_post_GSVA_6_thesis_grid.R",
  "7_post_GSVA_7_interaction.R",
  "7_post_GSVA_8_heatmap.R"
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
