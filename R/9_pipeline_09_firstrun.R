here::i_am("R/9_pipeline_09_firstrun.R")

########################
# Run all pipeline
########################

# Run scripts of pipeline
c(
  "9_pipeline_01_load.R",
  "9_pipeline_02_dge.R",
  "9_pipeline_03_plot_pre_DE.R",
  "9_pipeline_04_plot_post_DE.R",
  "9_pipeline_05_WGCNA.R",
  "9_pipeline_06_camera.R",
  "9_pipeline_07_GSVA.R",
  "9_pipeline_08_evaluate.R"
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
