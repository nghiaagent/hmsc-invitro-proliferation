here::i_am("R/6_post_camera_02_perform.R")

########################
# Perform camera on all contrasts
# Load function script beforehand
########################

# Import packages
library(DESeq2)
library(here)
library(tidyverse)

camera_all <- purrr::map(
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
)

# Save data
names(camera_all) <- colnames(contrasts)

saveRDS(
  camera_all,
  here::here(
    "output",
    "data_enrichment",
    "camera",
    "camera_all.RDS"
  )
)
