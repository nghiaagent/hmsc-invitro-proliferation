genenames_inhouse <- c(
  "ENO2",
  "NANOG",
  "NES",
  "TUBB3",
  "VIM",
  "ACAN",
  "HSPG2",
  "DCN",
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
  "HS3ST3A1",
  "HS3ST3B1",
  "NDST1",
  "NDST2",
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

geneIds_inhouse <- mapIds(org.Hs.eg.db,
                          keys = genenames_inhouse,
                          column = "ENTREZID",
                          keytype = "SYMBOL")

## Define genes in WNT signalling

msigdb_GOBP <- getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.entrez.gmt")

msigdb_GOMF <- getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.mf.v2023.2.Hs.entrez.gmt")

geneIds_Wnt <- c(msigdb_GOBP[["GOBP_CELL_CELL_SIGNALING_BY_WNT"]]@geneIds,
                 msigdb_GOBP[["GOBP_NON_CANONICAL_WNT_SIGNALING_PATHWAY"]]@geneIds,
                 msigdb_GOBP[["GOBP_CANONICAL_WNT_SIGNALING_PATHWAY"]]@geneIds) %>%
  unique()

names(geneIds_Wnt) <- mapIds(org.Hs.eg.db,
                             keys = geneIds_Wnt,
                             column = "SYMBOL",
                             keytype = "ENTREZID")

