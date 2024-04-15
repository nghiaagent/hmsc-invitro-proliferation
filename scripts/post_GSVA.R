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

quant_DGE_ENTREZ <- ExpressionSet(assayData = temp$E,
                                  phenoData = AnnotatedDataFrame(temp$targets),
                                  featureData = AnnotatedDataFrame(temp$genes))

