here::i_am("R/7_post_GSVA_6_thesis_grid.R")

########################
# Build plot grid for thesis
########################

# Import packages
library(conflicted)
library(cowplot)
library(GSVA)
library(here)
library(patchwork)
library(tidyverse)

# Load data
## Plots for passages
plots_poi_passage <- readRDS(
  here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "plots_poi_passage.RDS"
  )
) %>%
  unlist(recursive = FALSE)

# Create grid
## Get legend, reformat
plot_legend_passage <- cowplot::get_legend(
  plots_poi_passage[[1]] +
    ggplot2::guides(col = ggplot2::guide_legend(ncol = 3))
)

## Passages figure
plots_sel_poi_passage <- plots_poi_passage %>%
  purrr::map(\(plot) {
    plot +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 8)
      )
  })

plots_sel_poi_passage_horizontal <- patchwork::wrap_plots(
  plots_sel_poi_passage,
  byrow = FALSE
) /
  plot_legend_passage +
  patchwork::plot_layout(heights = c(20, 1))

plots_sel_poi_passage_vertical <- patchwork::wrap_plots(
  plots_sel_poi_passage,
  byrow = TRUE
) /
  plot_legend_passage +
  patchwork::plot_layout(heights = c(20, 1))

# Export plots
ggplot2::ggsave(
  filename = "pois_passage_horizontal.png",
  plot = plots_sel_poi_passage_horizontal,
  path = here::here(
    "output",
    "data_enrichment",
    "GSVA"
  ),
  scale = 0.8,
  width = 10,
  height = 11,
  units = "in",
  dpi = 144
)

ggplot2::ggsave(
  filename = "pois_passage_vertical.png",
  plot = plots_sel_poi_passage_vertical,
  path = here::here(
    "output",
    "data_enrichment",
    "GSVA"
  ),
  scale = 0.8,
  width = 10,
  height = 11,
  units = "in",
  dpi = 144
)
