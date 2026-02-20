here::i_am("R/9_pipeline_08_evaluate.R")

########################
# Run evaluation of methods
########################

# Run scripts of pipeline
c(
  "8_evaluate_batchcor.R",
  "8_evaluate_camera.R",
  "8_evaluate_clusterprofiler.R",
  "8_evaluate_glm.R",
  "8_plot_heatmap_evalbatchcor.R"
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
