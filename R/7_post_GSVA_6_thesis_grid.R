here::i_am("R/7_post_GSVA_6_thesis_grid.R")

########################
# Build plot grid for thesis
########################

# Import packages
library(here)
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
plot_legend_passage <- get_legend(
  plots_poi_passage[[1]] + guides(col = guide_legend(ncol = 3))
)

## Passages figure
plots_sel_poi_passage <- plots_poi_passage %>%
  map(\(plot) {
    plot +
      theme(
        legend.position = "none",
        plot.title = element_text(size = 8)
      )
  })

plots_sel_poi_passage_horizontal <- wrap_plots(
  plots_sel_poi_passage,
  byrow = FALSE
) /
  plot_legend_passage +
  plot_layout(heights = c(20, 1))

plots_sel_poi_passage_vertical <- wrap_plots(
  plots_sel_poi_passage,
  byrow = TRUE
) /
  plot_legend_passage +
  plot_layout(heights = c(20, 1))

# Export plots
ggsave(
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

ggsave(
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
