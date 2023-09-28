### Don't source this file by itself; call in from another file after running env_prep.R

source("./scripts/dge.R")


# Prepare dataset
# Swap ENSEMBL ID with ENTREZ ID
# Remove genes with no Entrez ID
# Remove genes with duplicated Entrez ID

quant_DGE_EGSEA <- quant_DGE_voom[!is.na(quant_DGE_voom$genes$ENTREZID),]
quant_DGE_EGSEA <- quant_DGE_EGSEA[!duplicated(quant_DGE_EGSEA$genes$ENTREZID),]
rownames(quant_DGE_EGSEA$E) <- quant_DGE_EGSEA$genes$ENTREZID

# Prepare gene sets

gs_annotations <- EGSEA::buildIdx(entrezIDs = quant_DGE_EGSEA$genes$ENTREZID,
                                  species = "human")

# Run EGSEA

egsea_passage <- egsea(
  voom.results = quant_DGE_EGSEA,
  contrasts = matrix_contrasts_passage,
  gs.annots = gs_annotations,
  symbolsMap = quant_DGE_EGSEA$genes[,c(2,3)],
  baseGSEAs = c("camera",
                "roast",
                "plage",
                "zscore",
                "gsva",
                "ssgsea",
                "ora",
                "fry"),
  report.dir = "./output/GSEA",
  sort.by = "avg.rank",
  report = TRUE
)
