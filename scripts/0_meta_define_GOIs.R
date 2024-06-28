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

geneIds_inhouse <- mapIds(org.Hs.eg.db,
                          keys = genenames_inhouse,
                          column = "ENTREZID",
                          keytype = "SYMBOL")

geneIds_rt2array <- mapIds(org.Hs.eg.db,
                           keys = genenames_rt2array,
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

## Define genes related to HSPGs (not already in in-house set)

msigdb_KEGG <- getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cp.kegg_legacy.v2023.2.Hs.entrez.gmt")

geneIds_HSPGs <- c(msigdb_KEGG[["KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_HEPARAN_SULFATE"]]@geneIds,
                   msigdb_GOBP[["GOBP_HEPARAN_SULFATE_PROTEOGLYCAN_BIOSYNTHETIC_PROCESS"]]@geneIds,
                   msigdb_GOMF[["GOMF_HEPARAN_SULFATE_PROTEOGLYCAN_BINDING"]]@geneIds) %>%
  unique() %>%
  setdiff(geneIds_inhouse)

names(geneIds_HSPGs) <- mapIds(org.Hs.eg.db,
                               keys = geneIds_HSPGs,
                               column = "SYMBOL",
                               keytype = "ENTREZID")
