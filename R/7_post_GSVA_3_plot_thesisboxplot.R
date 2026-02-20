here::i_am("R/7_post_GSVA_3_plot_thesisboxplot.R")

########################
# Load data
########################

# Import packages
library(DESeq2)
library(here)

fit_gsva <- readRDS(
  here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "GSVA_results.RDS"
  )
)

quant_gsva <- readRDS(
  here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "quant_GSVA.RDS"
  )
)

#
