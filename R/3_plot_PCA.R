# Draw PC1+2 PCA biplots for pre and post batch correction data
# Load data
quant_deseq2 <- readRDS(here::here(
  "output",
  "data_expression",
  "pre_DGE",
  "quant_cDNA_deseq.RDS"
))

rlog_deseq2 <- quant_deseq2 %>%
  rlog() %>%
  assay()

quant_deseq2_batchcor <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

rlog_deseq2_batchcor <- quant_deseq2_batchcor %>%
  rlog() %>%
  assay()

# Mutate quant objects to anonymise batch
# Create list of objects for plotting
list_rlog <- list(
  uncorrected = rlog_deseq2,
  corrected = rlog_deseq2_batchcor
)

# Get metadata table
quant_coldata <- colData(quant_deseq2)
quant_coldata$run_date <- quant_coldata$run_date %>%
  fct_anon()

# Define groupings for PCA plots
groupings <- c(
  "run_date",
  "cell_line",
  "timepoint_ID",
  "Treatment"
)

# Plot
## Get list of all PCA biplots
list_biplot_all <- map2(
  groupings,
  palette_pca,
  \(grouping, palette) {
    # Get list of PCA biplots
    list_plots <- map(
      list_rlog,
      \(rlog) {
        biplot <- rlog %>%
          pca(
            metadata = quant_coldata,
            removeVar = 0.9
          ) %>%
          PCAtools::biplot(
            colby = grouping,
            colkey = palette,
            legendPosition = "right",
            lab = NULL
          )
        # Return data
        return(biplot)
      }
    )

    # Organise into pairs
    grid_plots <- list_plots[[1]] +
      theme(legend.position = "none") +
      list_plots[[2]] +
      theme(legend.position = "none") +
      get_legend(list_plots[[1]]) +
      plot_layout(
        ncol = 3,
        widths = c(2, 2, 1)
      )
    # Return data
    return(grid_plots)
  }
)

plot_correlation <- map(
  list_rlog,
  \(rlog) {
    eigencorplot <- rlog %>%
      pca(
        metadata = quant_coldata,
        removeVar = 0.9
      ) %>%
      PCAtools::eigencorplot(
        metavars = c(
          "condition_ID",
          "timepoint_ID",
          "run_date",
          "Treatment",
          "Day",
          "Passage",
          "cell_line"
        ),
        col = hcl.colors(9, palette = "Blue-Red 2"),
        colCorval = "black"
      )
    # Return data
    return(eigencorplot)
  }
) %>%
  plot_grid(plotlist = .)

grid_biplot <- plot_grid(
  plotlist = list_biplot_all,
  ncol = 2,
  labels = "AUTO"
) %>%
  plot_grid(
    plot_correlation,
    labels = c("", "E"),
    ncol = 1,
    rel_heights = c(2, 1)
  )

# Selected plots (For paper)
## Batch corrected, colour = timepoint, shape = cell population
biplot_selected <- rlog_deseq2_batchcor %>%
  pca(metadata = quant_coldata) %>%
  PCAtools::biplot(
    colby = "timepoint_ID",
    colkey = c(
      P5D3 = palette_merge[[1]],
      P5D5 = palette_merge[[2]],
      P7D3 = palette_merge[[3]],
      P7D5 = palette_merge[[4]],
      P13D3 = palette_merge[[5]],
      P13D5 = palette_merge[[6]]
    ),
    colLegendTitle = "Timepoint",
    shape = "cell_line",
    shapeLegendTitle = "Cell population",
    legendPosition = "left",
    lab = NULL
  )

# Export plots
## All plots
ggsave(
  here::here(
    "output",
    "plots_pca",
    "biplots_all.png"
  ),
  grid_biplot,
  width = 16,
  height = 10,
  scale = 1.3
)

## Selected plots
ggsave(
  here::here(
    "output",
    "plots_pca",
    "biplots_selected.png"
  ),
  biplot_selected,
  width = 16,
  height = 8,
  scale = 0.5
)

# Random code: 3D PCA plot
list_pca <- list_rlog %>%
  map(\(x) {
    x %>%
      pca(metadata = quant_coldata, removeVar = 0.9)
  })

plot_ly(
  list_pca[["corrected"]]$rotated,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = list_pca[["corrected"]]$metadata$run_date
)
