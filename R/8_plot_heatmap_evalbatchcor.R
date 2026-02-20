here::i_am("R/8_plot_heatmap_evalbatchcor.R")

########################
# Plot significant DEGs as a heatmap
########################

# Import packages
library(ComplexHeatmap)
library(conflicted)
library(DESeq2)
library(here)
library(SummarizedExperiment)
library(tidyverse)

# Load data
quant_deseq2_nobatchcor <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2.RDS"
  )
)

quant_rlog_nobatchcor <- DESeq2::rlog(quant_deseq2_nobatchcor)

quant_deseq2_lrt <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_LRT.RDS"
  )
)

# Get significant genes
genes_significant <- quant_deseq2_lrt %>%
  DESeq2::results() %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  rownames()

quant_heatmap <- quant_rlog_nobatchcor[genes_significant, ]

rlog_heatmap <- quant_heatmap %>%
  SummarizedExperiment::assay() %>%
  t() %>%
  scale() %>%
  t()

# Set color scheme and breaks
## For gene expression (rlog z-scores)
col <- viridis::inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

# Create extra factors for splitting columns
split_cell_line <- SummarizedExperiment::colData(quant_heatmap)$cell_line

split_cell_line_passage <- stringr::str_c(
  SummarizedExperiment::colData(quant_heatmap)$cell_line,
  SummarizedExperiment::colData(quant_heatmap)$Passage,
  sep = "_"
) %>%
  factor(
    levels = c(
      "hMSC-20176_P5",
      "hMSC-20176_P7",
      "hMSC-20176_P13",
      "hMSC-21558_P5",
      "hMSC-21558_P7",
      "hMSC-21558_P13"
    )
  )

split_cell_line_day <- stringr::str_c(
  SummarizedExperiment::colData(quant_heatmap)$cell_line,
  SummarizedExperiment::colData(quant_heatmap)$Day,
  sep = "_"
) %>%
  factor(
    levels = c(
      "hMSC-20176_D3",
      "hMSC-20176_D5",
      "hMSC-21558_D3",
      "hMSC-21558_D5"
    )
  )

# Create vectors for ordering samples in specific cases of splits
order <- SummarizedExperiment::colData(quant_heatmap) %>%
  data.frame() %>%
  dplyr::arrange(cell_line, Passage, Day, Treatment) %>%
  rownames()

# Build annotation; include only necessary metadata
## Dataframe of annotation data
anno <- tibble::tibble(
  `Cell population` = SummarizedExperiment::colData(quant_heatmap)$cell_line,
  `Passage` = SummarizedExperiment::colData(quant_heatmap)$Passage,
  `Day` = SummarizedExperiment::colData(quant_heatmap)$Day,
  `Treatment` = SummarizedExperiment::colData(quant_heatmap)$Treatment %>%
    dplyr::case_match(
      "Untreated" ~ "Control",
      "Treated" ~ "Heparin"
    ),
  `Batch` = SummarizedExperiment::colData(quant_heatmap)$run_date %>%
    forcats::fct_recode()
)

## ComplexHeatmap metadata annotation object
anno_object <- ComplexHeatmap::HeatmapAnnotation(
  df = anno,
  which = "col",
  col = palette_heatmap,
  annotation_height = 0.6,
  annotation_legend_param = list(
    `Cell population` = list(
      nrow = 2,
      title = "Cell population",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    `Passage` = list(
      nrow = 3,
      title = "Passage",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    `Day` = list(
      nrow = 2,
      title = "Day",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    `Treatment` = list(
      nrow = 2,
      title = "Treatment",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    )
  )
)

# Create heatmap object
heatmap <- ComplexHeatmap::Heatmap(
  rlog_heatmap,
  name = "Gene\nZ-\nscore",
  col = circlize::colorRamp2(breaks, col),
  border = FALSE,

  # parameters for the colour-bar that represents gradient of expression
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    legend_width = unit(8, "cm"),
    legend_height = unit(5.0, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8, fontface = "bold")
  ),

  # row (gene) parameters
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  row_title = NULL,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 90,
  row_km = 3,
  row_km_repeats = 5,
  show_row_names = FALSE,

  # column (sample) parameters
  column_split = NULL,
  column_order = order,
  column_title = NULL,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_column_names = FALSE,

  # specify top and bottom annotations
  top_annotation = anno_object
)

# Export heatmap
png(
  file = "./output/plots_heatmap/normalised_expression_nobatchcor.png",
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

plot(heatmap)

dev.off()

# Redo the above with corrected data
# Plot significant DEGs as a heatmap
# Load data (rlog)
quant_deseq2_batchcor <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

quant_rlog_batchcor <- rlog(quant_deseq2_batchcor)

quant_deseq2_lrt <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_LRT.RDS"
  )
)

# Get significant genes
genes_significant <- quant_deseq2_lrt %>%
  results() %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.05) %>%
  rownames()

quant_heatmap <- quant_rlog_batchcor[genes_significant, ]

rlog_heatmap <- quant_heatmap %>%
  SummarizedExperiment::assay() %>%
  t() %>%
  scale() %>%
  t()

# Set color scheme and breaks
## For gene expression (rlog z-scores)
col <- viridis::inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

# Create extra factors for splitting columns
split_cell_line <- SummarizedExperiment::colData(quant_heatmap)$cell_line

split_cell_line_passage <- stringr::str_c(
  SummarizedExperiment::colData(quant_heatmap)$cell_line,
  SummarizedExperiment::colData(quant_heatmap)$Passage,
  sep = "_"
) %>%
  factor(
    levels = c(
      "hMSC-20176_P5",
      "hMSC-20176_P7",
      "hMSC-20176_P13",
      "hMSC-21558_P5",
      "hMSC-21558_P7",
      "hMSC-21558_P13"
    )
  )

split_cell_line_day <- stringr::str_c(
  SummarizedExperiment::colData(quant_heatmap)$cell_line,
  SummarizedExperiment::colData(quant_heatmap)$Day,
  sep = "_"
) %>%
  factor(
    levels = c(
      "hMSC-20176_D3",
      "hMSC-20176_D5",
      "hMSC-21558_D3",
      "hMSC-21558_D5"
    )
  )

# Create vectors for ordering samples in specific cases of splits
order <- SummarizedExperiment::colData(quant_heatmap) %>%
  data.frame() %>%
  dplyr::arrange(cell_line, Passage, Day, Treatment) %>%
  rownames()

# Build annotation; include only necessary metadata
## Dataframe of annotation data
anno <- tibble::tibble(
  `Cell population` = SummarizedExperiment::colData(quant_heatmap)$cell_line,
  `Passage` = SummarizedExperiment::colData(quant_heatmap)$Passage,
  `Day` = SummarizedExperiment::colData(quant_heatmap)$Day,
  `Treatment` = SummarizedExperiment::colData(quant_heatmap)$Treatment %>%
    dplyr::case_match(
      "Untreated" ~ "Control",
      "Treated" ~ "Heparin"
    ),
  `Batch` = SummarizedExperiment::colData(quant_heatmap)$run_date
)

## ComplexHeatmap metadata annotation object
anno_object <- ComplexHeatmap::HeatmapAnnotation(
  df = anno,
  which = "col",
  col = palette_heatmap,
  annotation_height = 0.6,
  annotation_legend_param = list(
    `Cell population` = list(
      nrow = 2,
      title = "Cell population",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    `Passage` = list(
      nrow = 3,
      title = "Passage",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    `Day` = list(
      nrow = 2,
      title = "Day",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    ),
    `Treatment` = list(
      nrow = 2,
      title = "Treatment",
      title_position = "topleft",
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 8, fontface = "bold")
    )
  )
)

# Create heatmap object
heatmap <- ComplexHeatmap::Heatmap(
  rlog_heatmap,
  name = "Gene\nZ-\nscore",
  col = circlize::colorRamp2(breaks, col),
  border = FALSE,

  # parameters for the colour-bar that represents gradient of expression
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    legend_width = unit(8, "cm"),
    legend_height = unit(5.0, "cm"),
    title_position = "topcenter",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8, fontface = "bold")
  ),

  # row (gene) parameters
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  row_title = NULL,
  row_title_side = "left",
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title_rot = 90,
  row_km = 3,
  row_km_repeats = 5,
  show_row_names = FALSE,

  # column (sample) parameters
  column_split = NULL,
  column_order = order,
  column_title = NULL,
  cluster_column_slices = FALSE,
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  show_column_names = FALSE,

  # specify top and bottom annotations
  top_annotation = anno_object
)

# Export heatmap
png(
  file = "./output/plots_heatmap/normalised_expression_batchcor.png",
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

plot(heatmap)

dev.off()
