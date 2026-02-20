here::i_am("R/3_plot_heatmap.R")

########################
# Plot significant DEGs as a heatmap
########################

# Import packages
library(DESeq2)
library(here)
library(SummarizedExperiment)
library(tidyverse)

# Load data (rlog and LRT data)
quant_rlog <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "rlog_deseq2.RDS"
  )
)

quant_deseq2_lrt <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_LRT.RDS"
  )
)

# Get significant genes
## Get list of significant genes
genes_significant <- quant_deseq2_lrt %>%
  results() %>%
  as.data.frame() %>%
  filter(padj < 0.05) %>%
  rownames()

## Filter heatmap to only significant genes
quant_heatmap <- quant_rlog[genes_significant, ]

## Get Z-scores for heatmap
rlog_heatmap <- quant_heatmap %>%
  assay() %>%
  t() %>%
  scale() %>%
  t()

# Set color scheme and breaks
## For gene expression (rlog z-scores)
col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

# Create extra factors for splitting columns
## The factors are:
## Cell population only
## Cell population + passage
## Cell population + day
split_cell_line <- colData(quant_heatmap)$cell_line

split_cell_line_passage <- str_c(
  colData(quant_heatmap)$cell_line,
  colData(quant_heatmap)$Passage,
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

split_cell_line_day <- str_c(
  colData(quant_heatmap)$cell_line,
  colData(quant_heatmap)$Day,
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
order <- quant_heatmap %>%
  colData() %>%
  data.frame() %>%
  arrange(cell_line, Passage, Day, Treatment) %>%
  rownames()

# Build annotation; include only necessary metadata
## Dataframe of annotation data
anno <- tibble(
  `Cell population` = colData(quant_heatmap)$cell_line,
  `Passage` = colData(quant_heatmap)$Passage,
  `Day` = colData(quant_heatmap)$Day,
  `Treatment` = colData(quant_heatmap)$Treatment %>%
    case_match(
      "Untreated" ~ "Control",
      "Treated" ~ "Heparin"
    )
)

## ComplexHeatmap metadata annotation object
anno_object <- HeatmapAnnotation(
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
heatmap <- Heatmap(
  rlog_heatmap,
  name = "Gene\nZ-\nscore",
  col = colorRamp2(breaks, col),
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
  column_split = split_cell_line_passage,
  column_order = order,
  column_title = NULL,
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,
  show_column_dend = FALSE,
  show_column_names = FALSE,

  # specify top and bottom annotations
  top_annotation = anno_object
)

# Export heatmap
png(
  file = "./output/plots_heatmap/normalised_expression_genes.png",
  width = 8,
  height = 12,
  units = "in",
  res = 600
)

plot(heatmap)

dev.off()
