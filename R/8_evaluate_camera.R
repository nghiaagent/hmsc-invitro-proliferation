here::i_am("R/8_evaluate_camera.R")

########################
# Source relevant scripts
########################

# Import packages
library(DESeq2)
library(here)
library(magrittr)
library(tidyverse)

source(here::here(
  "scripts",
  "6_post_camera_function.R"
))

# Run camera
camera_all_withigc <- purrr::map(
  seq_len(ncol(contrasts)),
  \(x) {
    run_camera(
      elist = rlog_camera,
      design = design,
      contrasts = contrasts,
      genesets = list_gmt,
      coef = x,
      inter.gene.cor = 0.01,
      sort = TRUE
    )
  }
) %>%
  set_names(colnames(contrasts))

camera_all_noigc <- purrr::map(
  seq_len(ncol(contrasts)),
  \(x) {
    run_camera(
      elist = rlog_camera,
      design = design,
      contrasts = contrasts,
      genesets = list_gmt,
      coef = x,
      inter.gene.cor = 0,
      sort = TRUE
    )
  }
) %>%
  set_names(colnames(contrasts))

# Plot
## Get relevant tables
list_camera <- camera_all_withigc %$%
  list(
    "P7vsP5_UT_D3" = .[["P7vsP5_UT_D3"]][["GOBP"]],
    "P13vsP7_UT_D3" = .[["P13vsP7_UT_D3"]][["GOBP"]],
    "P13vsP5_UT_D3" = .[["P13vsP5_UT_D3"]][["GOBP"]]
  ) %>%
  # Format data
  map(\(.results) {
    .results <- .results %>%
      dplyr::arrange(FDR) %>%
      dplyr::slice_head(n = 20) %>%
      rownames_to_column(var = "gs_name") %>%
      mutate(
        gs_name = gs_name %>%
          str_replace("GOBP_", "") %>%
          str_replace_all("_", " ") %>%
          factor(., levels = .)
      )

    # Return data
    return(.results)
  })

plots_camera <- list_camera %>%
  imap(\(.results, .name) {
    .plot <- .results %>%
      ggplot(aes(
        x = NGenes,
        y = gs_name,
        fill = FDR
      )) +
      geom_col() +
      scale_fill_continuous_sequential(palette = "Red-Blue") +
      theme_bw() +
      theme(
        axis.title.y = element_blank(),
        plot.title = element_text(face = "bold")
      ) +
      scale_y_discrete(
        labels = \(x) {
          str_wrap(
            x,
            width = 30,
            whitespace_only = TRUE
          )
        }
      ) +
      labs(
        x = "No. genes in set",
        title = .name,
      )
  })

grid_camera <- plots_camera %>%
  wrap_plots(nrow = 1)

# Save plot
ggsave(
  here::here(
    "output",
    "plots_QC",
    "camera eval.png"
  ),
  grid_camera,
  width = 16,
  height = 9,
  scale = 1.1
)
