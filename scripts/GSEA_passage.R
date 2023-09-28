### Don't source this file by itself; call in from another file after running env_prep.R

source("./scripts/dge.R")

# Prepare gene sets

gs_annotations <- EGSEA::buildIdx(entrezIDs = quant_DGE_voom$genes$ENTREZID,
                species = "human")

# Prepare dataset
# Swap ENSEMBL ID with ENTREZ ID

rownames(quant_DGE_voom$E) <- quant_DGE_voom$genes$ENTREZID

# Run EGSEA

egsea_passage <- egsea(
  voom.results = quant_DGE_voom,
  contrasts = matrix_contrasts_passage,
  gs.annots = gs_annotations,
  symbolsMap = quant_DGE_voom$genes[,c(2,3)],
  report.dir = "./outputs/GSEA",
  sort.by = "avg.rank",
  report = TRUE
)
