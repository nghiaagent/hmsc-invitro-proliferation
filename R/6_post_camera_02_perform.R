here::i_am("R/6_post_camera_02_perform.R")

########################
# Perform camera on all contrasts
########################

# Import packages
library(conflicted)
library(DESeq2)
library(here)
library(limma)
library(tidyverse)

# Make sure to source the prep script first
# source(here::here("R/6_post_camera_01_prepare.R"))

# Perform camera on all contrasts
camera_all <- purrr::map(
  seq_len(ncol(contrasts)),
  \(x) {
    run_camera(
      elist = rlog_camera,
      design = design,
      contrasts = contrasts,
      genesets = list_gmt_camera,
      coef = x,
      inter.gene.cor = 0.01,
      sort = TRUE
    )
  }
)

names(camera_all) <- colnames(contrasts)

# Save data
saveRDS(
  camera_all,
  here::here(
    "output",
    "data_enrichment",
    "camera",
    "camera_all.RDS"
  )
)
