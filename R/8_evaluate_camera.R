here::i_am("R/8_evaluate_camera.R")

########################
# Evaluate CAMERA
########################

# Import packages
library(conflicted)
library(DESeq2)
library(here)
library(magrittr)
library(patchwork)
library(tidyverse)

# Source relevant scripts
source(here::here(
  "R",
  "6_post_camera_01_prepare.R"
))

# Run camera
## With inter-gene correlation
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
  magrittr::set_names(colnames(contrasts))

## Without inter-gene correlation
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
  magrittr::set_names(colnames(contrasts))

# Plot
## Get relevant tables
list_camera <- camera_all_withigc %$%
  list(
    "P7vsP5_UT_D3" = .[["P7vsP5_UT_D3"]][["GOBP"]],
    "P13vsP7_UT_D3" = .[["P13vsP7_UT_D3"]][["GOBP"]],
    "P13vsP5_UT_D3" = .[["P13vsP5_UT_D3"]][["GOBP"]]
  ) %>%
  # Format data
  purrr::map(\(.results) {
    .results <- .results %>%
      dplyr::arrange(FDR) %>%
      dplyr::slice_head(n = 20) %>%
      tibble::rownames_to_column(var = "gs_name") %>%
      dplyr::mutate(
        gs_name = gs_name %>%
          stringr::str_replace("GOBP_", "") %>%
          stringr::str_replace_all("_", " ") %>%
          factor(., levels = .)
      )

    # Return data
    return(.results)
  })

## Build plots
plots_camera <- list_camera %>%
  purrr::imap(\(.results, .name) {
    .plot <- .results %>%
      ggplot2::ggplot(ggplot2::aes(
        x = NGenes,
        y = gs_name,
        fill = FDR
      )) +
      ggplot2::geom_col() +
      colorspace::scale_fill_continuous_sequential(palette = "Red-Blue") +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::scale_y_discrete(
        labels = \(x) {
          stringr::str_wrap(
            x,
            width = 30,
            whitespace_only = TRUE
          )
        }
      ) +
      ggplot2::labs(
        x = "No. genes in set",
        title = .name,
      )
  })

grid_camera <- plots_camera %>%
  patchwork::wrap_plots(nrow = 1)

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
