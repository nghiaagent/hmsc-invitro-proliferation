here::i_am("R/9_pipeline_06_camera.R")

########################
# Run CAMERA gene set analysis
########################

# Run scripts of pipeline
c(
  "6_post_camera_01_prepare.R",
  "6_post_camera_02_perform.R",
  "6_post_camera_03_plot.R",
  "6_post_camera_04_plot_thesis_grid.R"
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
