here::i_am("R/0_define_GOIs.R")

########################
# Define in-house GOIs
########################

# Import packages
library(GSEABase)
library(here)
library(org.Hs.eg.db)
library(tidyverse)

# Define gene names (non-HSPGs and HSPGs)
genenames_inhouse <- c(
  "ENO2",
  "NANOG",
  "NES",
  "TUBB3",
  "VIM",
  "GPC1",
  "GPC4",
  "SDC1",
  "SDC2",
  "SDC3",
  "SDC4",
  "EXT1",
  "EXT2",
  "GLCE",
  "HS2ST1",
  "HS3ST1",
  "HS3ST2",
  "HS3ST3A1",
  "HS3ST3B1",
  "HS3ST4",
  "HS3ST6",
  "HS6ST1",
  "HS6ST2",
  "HS6ST3",
  "HPSE2",
  "NDST1",
  "NDST2",
  "NDST3",
  "NDST4",
  "ACTA2",
  "CD44",
  "COL1A1",
  "COL1A2",
  "RUNX2",
  "CHST11",
  "SOX1",
  "MAP2",
  "S100B",
  "POU5F1"
)

genenames_inhouse_hspgs <- c(
  "VIM",
  "GPC1",
  "GPC4",
  "SDC1",
  "SDC2",
  "SDC3",
  "SDC4",
  "EXT1",
  "EXT2",
  "GLCE",
  "HS2ST1",
  "HS3ST1",
  "HS3ST2",
  "HS3ST3A1",
  "HS3ST3B1",
  "HS3ST4",
  "HS3ST6",
  "HS6ST1",
  "HS6ST2",
  "HS6ST3",
  "HPSE2",
  "NDST1",
  "NDST2",
  "NDST3",
  "NDST4"
)

## Map gene names to ENTREZIDs
geneids_inhouse <- mapIds(
  org.Hs.eg.db,
  keys = genenames_inhouse,
  column = "ENTREZID",
  keytype = "SYMBOL"
)

geneids_inhouse_hspgs <- mapIds(
  org.Hs.eg.db,
  keys = genenames_inhouse_hspgs,
  column = "ENTREZID",
  keytype = "SYMBOL"
)

# Define AD-related genes from QIAgen panel (for SIP)
genenames_rt2array <- c(
  "A2M",
  "ABCA1",
  "ACHE",
  "ADAM10",
  "ADAM9",
  "APBA1",
  "APBA3",
  "APBB1",
  "APBB2",
  "APH1A",
  "APLP1",
  "APLP2",
  "APOA1",
  "APOE",
  "APP",
  "BACE1",
  "BACE2",
  "BCHE",
  "BDNF",
  "CAPN1",
  "CASP3",
  "CASP4",
  "CDK1",
  "CDK5",
  "CDKL1",
  "CHAT",
  "CLU",
  "CTSB",
  "CTSC",
  "CTSD",
  "CTSG",
  "CTSL",
  "EP300",
  "ERN1",
  "GAP43",
  "GNAO1",
  "GNAZ",
  "GNB1",
  "GNB2",
  "GNB4",
  "GNB5",
  "GNG11",
  "GNG3",
  "GNG4"
)

## Map gene names to ENTREZIDs
geneids_rt2array <- mapIds(
  org.Hs.eg.db,
  keys = genenames_rt2array,
  column = "ENTREZID",
  keytype = "SYMBOL"
)

# Load external GMTs
msigdb_gobp <- getGmt(
  con = here::here(
    "input",
    "genesets",
    "msigdb_v2023.2.Hs_GMTs",
    "c5.go.bp.v2023.2.Hs.entrez.gmt"
  )
)

msigdb_gomf <- getGmt(
  con = here::here(
    "input",
    "genesets",
    "msigdb_v2023.2.Hs_GMTs",
    "c5.go.mf.v2023.2.Hs.entrez.gmt"
  )
)

msigdb_kegg <- getGmt(
  con = here::here(
    "input",
    "genesets",
    "msigdb_v2023.2.Hs_GMTs",
    "c2.cp.kegg_legacy.v2023.2.Hs.entrez.gmt"
  )
)

list_gmt <- list(
  "h" = "msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.entrez.gmt",
  "c2_cgp" = "msigdb_v2023.2.Hs_GMTs/c2.cgp.v2023.2.Hs.entrez.gmt",
  "c2_cp" = "msigdb_v2023.2.Hs_GMTs/c2.cp.v2023.2.Hs.entrez.gmt",
  "GOBP" = "msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.entrez.gmt",
  "GOCC" = "msigdb_v2023.2.Hs_GMTs/c5.go.cc.v2023.2.Hs.entrez.gmt",
  "GOMF" = "msigdb_v2023.2.Hs_GMTs/c5.go.mf.v2023.2.Hs.entrez.gmt",
  "WGCNA" = "gcn_sets_WGCNA.gmt"
) %>%
  map(\(filename) {
    here::here(
      "input",
      "genesets",
      filename
    )
  }) %>%
  map(\(x) getGmt(con = x))

## Define genes in WNT signalling
geneids_wnt <- c(
  msigdb_gobp[["GOBP_CELL_CELL_SIGNALING_BY_WNT"]]@geneIds,
  msigdb_gobp[["GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY"]]@geneIds,
  msigdb_gobp[["GOBP_CANONICAL_WNT_SIGNALING_PATHWAY"]]@geneIds
)

names(geneids_wnt) <- mapIds(
  org.Hs.eg.db,
  keys = geneids_wnt,
  column = "SYMBOL",
  keytype = "ENTREZID"
)

## Define genes related to HSPGs (not already in in-house set)
geneids_hspgs <- c(
  msigdb_kegg[["KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_HEPARAN_SULFATE"]]@geneIds,
  msigdb_gobp[[
    "GOBP_HEPARAN_SULFATE_PROTEOGLYCAN_BIOSYNTHETIC_PROCESS"
  ]]@geneIds,
  msigdb_gomf[["GOMF_HEPARAN_SULFATE_PROTEOGLYCAN_BINDING"]]@geneIds
)

names(geneids_hspgs) <- mapIds(
  org.Hs.eg.db,
  keys = geneids_hspgs,
  column = "SYMBOL",
  keytype = "ENTREZID"
)

# Combine into one lists
geneids_goi <- c(
  geneids_inhouse,
  geneids_hspgs,
  geneids_inhouse_hspgs,
  geneids_rt2array,
  geneids_wnt
) %>%
  unique()

geneids_goi_limited <- c(
  geneids_inhouse,
  geneids_inhouse_hspgs,
  geneids_hspgs
) %>%
  unique()

# Define genes related to NF-kB signalling
geneids_nfkb <- c(
  list_gmt[["h"]][["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]@geneIds
) %>%
  unique() %>%
  set_names(mapIds(
    org.Hs.eg.db,
    keys = .,
    column = "SYMBOL",
    keytype = "ENTREZID"
  ))

# Define genes related to TGF-b signallings
geneids_tgfb <- c(
  list_gmt[["h"]][["HALLMARK_TGF_BETA_SIGNALING"]]@geneIds
) %>%
  unique() %>%
  set_names(mapIds(
    org.Hs.eg.db,
    keys = .,
    column = "SYMBOL",
    keytype = "ENTREZID"
  ))
