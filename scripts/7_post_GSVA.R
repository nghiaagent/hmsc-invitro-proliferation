# Run relevant scripts beforehand

# env_prep
# dge_cellpops_as_fixed

# Convert quant_DGE_voom to a compatible format
## Convert EList to ExpressionSet
## Convert ENSEMBL ID to ENTREZ ID

temp <- quant_DGE_voom

rownames(temp$genes) <- temp$genes$GENEID

order <- order(fit$Amean, decreasing = TRUE)

temp <- temp[order,]

temp <- temp[complete.cases(temp$genes),]

temp <- temp[!duplicated(temp$genes$ENTREZID),]

rownames(temp$E) <- temp$genes$ENTREZID

rownames(temp$genes) <- temp$genes$ENTREZID

temp <- temp[order(temp$genes$GENEID, decreasing = FALSE),]

temp$genes <- temp$genes %>% select(!GENEID)

quant_DGE_ESet <- ExpressionSet(
  assayData = temp$E,
  phenoData = AnnotatedDataFrame(temp$targets),
  featureData = AnnotatedDataFrame(temp$genes)
)

rm(temp)

# Build ontology gene sets

## GO

### GOBP
### Represented by MSigDB c5/GO/BP

msigdb_GOBP <-
  getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.entrez.gmt")

### GOMF
### Represented by MSigDB c5/GO/MF

msigdb_GOMF <-
  getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.mf.v2023.2.Hs.entrez.gmt")

### GOCC
### Represented by MSigDB c5/GO/CC

msigdb_GOCC <-
  getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.cc.v2023.2.Hs.entrez.gmt")

## MSigDB

msigdb_h  <-
  getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.entrez.gmt")

msigdb_c2 <-
  getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.all.v2023.2.Hs.entrez.gmt")

msigdb_c3 <-
  getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c3.all.v2023.2.Hs.entrez.gmt")

## ReactomePA
## Represented by MSigDB c2/CP/Reactome

msigdb_reactome <-
  getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cp.reactome.v2023.2.Hs.entrez.gmt")

## KEGG
## Represented by MSigDB c2/CP/KEGG_LEGACY

msigdb_KEGG <-
  getGmt(con = "./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cp.kegg_legacy.v2023.2.Hs.entrez.gmt")

## GCN gene sets

GCN_sets <- getGmt(con = file.path(".",
                                   "input",
                                   "genesets",
                                   "GCN_sets.gmt"))

# Perform GSVA

## GO

fit_GSVA_GOBP <- gsvaParam(quant_DGE_ESet,
                           msigdb_GOBP,
                           minSize = 5,
                           maxSize = 500) %>%
  gsva() %>%
  lmFit(., design) %>%
  eBayes() %>%
  contrasts.fit(., matrix_contrasts) %>%
  eBayes()

fit_GSVA_GOMF <- gsvaParam(quant_DGE_ESet,
                           msigdb_GOMF,
                           minSize = 5,
                           maxSize = 500) %>%
  gsva() %>%
  lmFit(., design) %>%
  eBayes() %>%
  contrasts.fit(., matrix_contrasts) %>%
  eBayes()


fit_GSVA_GOCC <- gsvaParam(quant_DGE_ESet,
                           msigdb_GOCC,
                           minSize = 5,
                           maxSize = 500) %>%
  gsva() %>%
  lmFit(., design) %>%
  eBayes() %>%
  contrasts.fit(., matrix_contrasts) %>%
  eBayes()

## MSigDB

fit_GSVA_h <- gsvaParam(quant_DGE_ESet,
                        msigdb_h,
                        minSize = 5,
                        maxSize = 500) %>%
  gsva() %>%
  lmFit(., design) %>%
  eBayes() %>%
  contrasts.fit(., matrix_contrasts) %>%
  eBayes()

fit_GSVA_c2 <- gsvaParam(quant_DGE_ESet,
                         msigdb_c2,
                         minSize = 5,
                         maxSize = 500) %>%
  gsva() %>%
  lmFit(., design) %>%
  eBayes() %>%
  contrasts.fit(., matrix_contrasts) %>%
  eBayes()

fit_GSVA_c3 <- gsvaParam(quant_DGE_ESet,
                         msigdb_c3,
                         minSize = 5,
                         maxSize = 500) %>%
  gsva() %>%
  lmFit(., design) %>%
  eBayes() %>%
  contrasts.fit(., matrix_contrasts) %>%
  eBayes()

## Reactome

fit_GSVA_reactome <- gsvaParam(quant_DGE_ESet,
                               msigdb_reactome,
                               minSize = 5,
                               maxSize = 500) %>%
  gsva() %>%
  lmFit(., design) %>%
  eBayes() %>%
  contrasts.fit(., matrix_contrasts) %>%
  eBayes()

## KEGG

fit_GSVA_KEGG <- gsvaParam(quant_DGE_ESet,
                           msigdb_KEGG,
                           minSize = 5,
                           maxSize = 500) %>%
  gsva() %>%
  lmFit(., design) %>%
  eBayes() %>%
  contrasts.fit(., matrix_contrasts) %>%
  eBayes()

fit_GSVA_GCN <- gsvaParam(quant_DGE_ESet,
                          GCN_sets,
                          minSize = 5,
                          maxSize = Inf) %>%
  gsva() %>%
  lmFit(., design) %>%
  eBayes() %>%
  contrasts.fit(., matrix_contrasts) %>%
  eBayes()

# Combine fit objects into list

fit_GSVA <- list(
  "GOBP" = fit_GSVA_GOBP,
  "GOCC" = fit_GSVA_GOCC,
  "GOMF" = fit_GSVA_GOMF,
  "h" = fit_GSVA_h,
  "c2" = fit_GSVA_c2,
  "c3" = fit_GSVA_c3,
  "KEGG" = fit_GSVA_KEGG,
  "Reactome" = fit_GSVA_reactome,
  "GCN" = fit_GSVA_GCN
)

# Make GSVA NES object for plotting

## GO

quant_GSVA_GOBP <- gsvaParam(quant_DGE_ESet,
                             msigdb_GOBP,
                             minSize = 5,
                             maxSize = 500) %>%
  gsva()

quant_GSVA_GOMF <- gsvaParam(quant_DGE_ESet,
                             msigdb_GOMF,
                             minSize = 5,
                             maxSize = 500) %>%
  gsva()


quant_GSVA_GOCC <- gsvaParam(quant_DGE_ESet,
                             msigdb_GOCC,
                             minSize = 5,
                             maxSize = 500) %>%              
  gsva()


quant_GSVA_h <- gsvaParam(quant_DGE_ESet,
                             msigdb_h,
                             minSize = 5,
                             maxSize = 500) %>%              
  gsva()


quant_GSVA_c2 <- gsvaParam(quant_DGE_ESet,
                           msigdb_c2,
                           minSize = 5,
                           maxSize = 500) %>%
  gsva()

quant_GSVA_c3 <- gsvaParam(quant_DGE_ESet,
                           msigdb_c3,
                           minSize = 5,
                           maxSize = 500) %>%
  gsva()

## Reactome

quant_GSVA_reactome <- gsvaParam(quant_DGE_ESet,
                                 msigdb_reactome,
                                 minSize = 5,
                                 maxSize = 500) %>%
  gsva()

## KEGG

quant_GSVA_KEGG <- gsvaParam(quant_DGE_ESet,
                             msigdb_KEGG,
                             minSize = 5,
                             maxSize = 500) %>%
  gsva()

## GCN

quant_GSVA_GCN <- gsvaParam(quant_DGE_ESet,
                             GCN_sets,
                             minSize = 5,
                             maxSize = Inf) %>%
  gsva()


# Combine fit objects into list

quant_GSVA <- list(
  "GOBP" = quant_GSVA_GOBP,
  "GOCC" = quant_GSVA_GOCC,
  "GOMF" = quant_GSVA_GOMF,
  "h" = quant_GSVA_h,
  "c2" = quant_GSVA_c2,
  "c3" = quant_GSVA_c3,
  "KEGG" = quant_GSVA_KEGG,
  "Reactome" = quant_GSVA_reactome,
  "GCN" = quant_GSVA_GCN
)

# Save data

saveRDS(fit_GSVA,
        "./output/data_enrichment/GSVA/GSVA_results.RDS")

saveRDS(quant_GSVA,
        "./output/data_enrichment/GSVA/quant_GSVA.RDS")

# Summarise DEG list

## Interaction terms (Treat x Timepoint)

lapply(fit_GSVA, decideTests) %>%
  lapply(summary) %>%
  lapply(function(x)
    x[, c(25:30, 37:39)])

## Venn diagrams (Treat x Day)

lapply(fit_GSVA, decideTests) %>%
  lapply(summary) %>%
  lapply(function(x)
    x[, c(7:12)])

## Venn diagrams (Treat x Passage)

lapply(fit_GSVA, decideTests) %>%
  lapply(summary) %>%
  lapply(function(x)
    x[, c(13:24)])

#
# # Target stats analysis to desired gene sets
#
# ## Get names of desired gene sets
#
# sets_interest <- c(
#   msigdb_c2[grepl("_SULFATE_", names(msigdb_c2))] %>% names(),
#   msigdb_GOBP[grepl("_SULFATE_", names(msigdb_GOBP))] %>% names(),
#   msigdb_GOCC[grepl("_SULFATE_", names(msigdb_GOCC))] %>% names(),
#   msigdb_GOMF[grepl("_SULFATE_", names(msigdb_GOMF))] %>% names()
# )
#
# ## Interaction terms (Treat x Timepoint)
#
# lapply(fit_GSVA_all,
#        function(x) {
#          x[rownames(x$coefficients) %in% sets_interest,]
#        }) %>%
#   lapply(decideTests) %>%
#   lapply(summary) %>%
#   lapply(function(x) {
#     x[, c(25:30, 37:39)]
#   })