here::i_am("R/x_dge_GOIs.R")

########################
# Load data
########################

# Import packages
library(here)
library(tidyverse)

source("./scripts/dge_cellpops_as_fixed.R")

# Subset for genes based on my stage 2.

list_KEGG_pathways <- c(
  "hsa00534",
  "hsa04110",
  "hsa04010",
  "hsa04310",
  "hsa04512",
  "hsa04514",
  "hsa04550",
  "hsa05205"
)

genes_inhouse <- c(
  "ENO2",
  "NANOG",
  "NES",
  "TUBB3",
  "VIM",
  "ACAN",
  "HSPG2",
  "DCN",
  "GPC1",
  "GPC2",
  "GPC3",
  "GPC4",
  "GPC5",
  "GPC6",
  "SDC1",
  "SDC2",
  "SDC3",
  "SDC4",
  "EXT1",
  "EXT2",
  "GLCE",
  "HS2ST1",
  "HS3ST3A1",
  "HS3ST3B1",
  "EXT1",
  "EXT2",
  "GLCE",
  "NDST1",
  "NDST2",
  "ACTA2",
  "CD44",
  "COL1A1",
  "COL1A2",
  "RUNX2",
  "CHST11",
  "SOX1"
)

genes_E2F <- c(
  "E2F1",
  "E2F2",
  "E2F3",
  "E2F4",
  "E2F5",
  "E2F6",
  "E2F7",
  "E2F8"
)
#
# genes_kegg <- getGeneKEGGLinks(species="hsa")
#
# genes_kegg <- filter(genes_kegg,
#                      PathwayID %in% list_KEGG_pathways)
#
# genes_fit_small <- intersect(fit_contrasts$genes$ENTREZID,
#                              genes_kegg$GeneID)

entrezid_inhouse <- fit_contrasts$genes[
  fit_contrasts$genes$GENENAME %in% genes_E2F,
]$ENTREZID
#
# entrezid_kegg <- fit_contrasts$genes[fit_contrasts$genes$ENTREZID %in% genes_fit_small,]$ENTREZID

entrezid_small <- entrezid_inhouse

# Filter fit_contrasts for GOIs. There are 701 GOIs.

fit_small <- fit_contrasts[fit_contrasts$genes$ENTREZID %in% entrezid_small, ]
