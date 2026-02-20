here::i_am("R/4_dge_plot_GOI_boxplot_thesis_grid.R")

########################
# Plot GOIs in a grid for thesis using cowplot
########################

# Import packages
library(here)
library(tidyverse)

# Load data
## Plots for passages
plots_goi_passage <- readRDS(
  file = here::here(
    "output",
    "plots_boxplot_counts",
    "plots_goi_passage.RDS"
  )
)

## Plots for treatments
plots_goi_treat <- readRDS(
  file = here::here(
    "output",
    "plots_boxplot_counts",
    "plots_goi_treat.RDS"
  )
)

# Define genes for plots
genes_sel_hspgs <- c(
  "SDC1",
  "SDC2",
  "SDC3",
  "SDC4",
  "GPC1",
  "GPC4",
  "NDST1",
  "NDST2",
  "GLCE",
  "HS2ST1",
  "HS3ST3A1",
  "HS3ST3B1",
  "EXT1",
  "EXTL1",
  "EXTL2",
  "EXTL3"
)

genes_sel_markers <- c(
  "NES",
  "SOX1",
  "TUBB3",
  "ENO2",
  "ACTA2",
  "COL1A1",
  "COL1A2",
  "RUNX2",
  "CD44"
)

# Create grid
## Get legend, reformat
plot_legend_passage <- get_legend(
  plots_goi_passage[[1]] +
    guides(col = guide_legend(ncol = 3))
)

plot_legend_treat <- get_legend(
  plots_goi_treat[[1]] +
    guides(col = guide_legend(ncol = 3)) +
    theme(legend.title = element_blank())
)

## Passages HSPG figure
plots_sel_passage_hspgs <- plots_goi_passage[genes_sel_hspgs] %>%
  map(\(plot) plot <- plot + theme(legend.position = "none"))

plots_sel_passage_hspgs <- wrap_plots(plots_sel_passage_hspgs) /
  plot_legend_passage +
  plot_layout(heights = c(20, 1))

## Treatment HSPG figure
plots_sel_treat_hspgs <- plots_goi_treat[genes_sel_hspgs] %>%
  map(\(plot) plot <- plot + theme(legend.position = "none"))

plots_sel_treat_hspgs <- wrap_plots(plots_sel_treat_hspgs) /
  plot_legend_treat +
  plot_layout(heights = c(20, 1))

# Export plots
ggsave(
  filename = str_c("gois_passage_hspgs.png", sep = ""),
  plot = plots_sel_passage_hspgs,
  path = here::here(
    ".",
    "output",
    "plots_boxplot_counts"
  ),
  scale = 1,
  width = 10,
  height = 11,
  units = "in",
  dpi = 144
)

ggsave(
  filename = str_c("gois_treat_hspgs.png", sep = ""),
  plot = plots_sel_treat_hspgs,
  path = here::here(
    ".",
    "output",
    "plots_boxplot_counts"
  ),
  scale = 1,
  width = 10,
  height = 11,
  units = "in",
  dpi = 144
)

## Passages hMSC figure
plots_sel_passage_markers <- plots_goi_passage[genes_sel_markers] %>%
  map(\(plot) plot <- plot + theme(legend.position = "none"))

plots_sel_passage_markers <- wrap_plots(plots_sel_passage_markers) /
  plot_legend_passage +
  plot_layout(heights = c(20, 1))

## Treatment hMSC figure
plots_sel_treat_markers <- plots_goi_treat[genes_sel_markers] %>%
  map(\(plot) plot <- plot + theme(legend.position = "none"))

plots_sel_treat_markers <- wrap_plots(plots_sel_treat_markers) /
  plot_legend_treat +
  plot_layout(heights = c(20, 1))

# Export plots
ggsave(
  filename = str_c("gois_passage_markers.png", sep = ""),
  plot = plots_sel_passage_markers,
  path = here::here(
    ".",
    "output",
    "plots_boxplot_counts"
  ),
  scale = 0.8,
  width = 10,
  height = 11,
  units = "in",
  dpi = 144
)

ggsave(
  filename = str_c("gois_treat_markers.png", sep = ""),
  plot = plots_sel_treat_markers,
  path = here::here(
    ".",
    "output",
    "plots_boxplot_counts"
  ),
  scale = 0.8,
  width = 10,
  height = 11,
  units = "in",
  dpi = 144
)
