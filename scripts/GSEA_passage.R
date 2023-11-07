### Don't source this file by itself; call in from another file after running env_prep.R

source("./scripts/dge_cellpops_as_fixed.R")


# Prepare dataset
# Swap ENSEMBL ID with ENTREZ ID
# Remove genes with no Entrez ID
# Remove genes with duplicated Entrez ID

quant_DGE_EGSEA <- quant_DGE_voom[!is.na(quant_DGE_voom$genes$ENTREZID),]
quant_DGE_EGSEA <- quant_DGE_EGSEA[!duplicated(quant_DGE_EGSEA$genes$ENTREZID),]
rownames(quant_DGE_EGSEA$E) <- quant_DGE_EGSEA$genes$ENTREZID

# Prepare gene sets

gs_annotations <- EGSEA::buildIdx(entrezIDs = quant_DGE_EGSEA$genes$ENTREZID,
                                  species = "human",
                                  msigdb.gsets = c("h",
                                                   "c2",
                                                   "c5"))

## Remove problematic KEGG pathways

problematic_pathways <- c("Retrograde endocannabinoid signaling",
                          "Apoptosis - multiple species",
                          "MicroRNAs in cancer",
                          "Mucin type O-Glycan biosynthesis",
                          "Glycosphingolipid biosynthesis - lacto and neolacto series",
                          "Other glycan degradation",
                          "Biosynthesis of amino acids",
                          "Other types of O-glycan biosynthesis",
                          "Carbon metabolism",
                          "2-Oxocarboxylic acid metabolism",
                          "Biosynthesis of unsaturated fatty acids")

sel <- which(names(gs_annotations[["kegg"]]@original) %in% problematic_pathways)

gs_annotations[["kegg"]]@original <- gs_annotations[["kegg"]]@original[-sel] 
gs_annotations[["kegg"]]@idx <- gs_annotations[["kegg"]]@idx[-sel] 
gs_annotations[["kegg"]]@anno <- gs_annotations[["kegg"]]@anno[-sel, ]

# Run EGSEA

egsea_treatment <- egsea(
  voom.results = quant_DGE_EGSEA,
  contrasts = matrix_contrasts,
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
  sort.by = "avg.rank",
  report.dir = "./egsea_treatment",
  report = TRUE
)
