here::i_am("R/6_post_camera_03_plot.R")

########################
# Build plots for camera
########################

# Import packages
library(conflicted)
library(DESeq2)
library(here)
library(tidyverse)

# Load data
results <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2_lfcshrink.RDS"
  )
)

list_gmt_camera <- readRDS(
  file = here::here("output/data_enrichment/camera/list_gmt_camera.rds")
)

# Plot barcode plots for Hallmark gene sets
## Plot multiple plots using a nested map
## Inner map: Create barcodeplot of all gene sets on provided contrast
## Outer map: Perform inner map on all contrasts
list_barcodeplots <- purrr::imap(
  results,
  \(results, name_contrast) {
    # Set and create output directory
    path_output <- here::here(
      "output",
      "data_enrichment",
      "camera",
      name_contrast,
      "barcodeplots"
    )

    if (!dir.exists(path_output)) {
      dir.create(path_output, recursive = TRUE)
    }

    # Create and draw plot
    purrr::imap(
      list_gmt_camera$h,
      \(genesets, name_geneset) {
        # Set up output device
        png(
          filename = here::here(
            path_output,
            stringr::str_c(name_geneset, ".png", sep = "")
          ),
          res = 288,
          width = 4,
          height = 3,
          units = "in"
        )

        par(mar = c(4.5, 3, 1.5, 1) - 0.5)

        # Draw plot
        limma::barcodeplot(
          results$log2FoldChange,
          index = genesets,
          alpha = 0.2,
          labels = c("", ""),
          xlab = "Log2 fold change",
          cex.lab = 0.6,
          cex.axis = 0.6,
          cex.main = 0.6,
          cex.sub = 0.6
        )

        # Export plot and object
        x <- recordPlot()
        dev.off()
        return(x)
      },
      .progress = TRUE
    )
  },
  .progress = TRUE
)

# Save data
saveRDS(
  list_barcodeplots,
  file = here::here(
    "output",
    "data_enrichment",
    "camera",
    "camera_barcodeplots.RDS"
  )
)
