# Load data

source("./scripts/dge_cellpops_as_fixed.R")

# Prepare dataset
# Swap ENSEMBL ID with ENTREZ ID
# Remove genes with no Entrez ID
# Remove genes with duplicated Entrez ID
# Re-run DGE with only 3 passage contrasts

quant_DGE_EGSEA <- quant_DGE_voom[!is.na(quant_DGE_voom$genes$ENTREZID),]
quant_DGE_EGSEA <- quant_DGE_EGSEA[!duplicated(quant_DGE_EGSEA$genes$ENTREZID),]
rownames(quant_DGE_EGSEA$E) <- quant_DGE_EGSEA$genes$ENTREZID

## Apply limma model fit

fit <- lmFit(quant_DGE_EGSEA,
             design) %>%
  eBayes()

# Define posthoc tests
## Between treatments, at each passage
## Each coefficient calculates the difference between that condition and the baseline (P5D3Untreated)

matrix_contrasts <- makeContrasts(
  P13vsP5 = condition_IDP13D3Untreated - 0,
  P13vsP7 = condition_IDP13D3Untreated - condition_IDP7D3Untreated,
  levels = design
)

# Test for desired contrasts

fit_contrasts <- contrasts.fit(fit,
                               matrix_contrasts) %>%
  eBayes()

summary(decideTests(fit_contrasts))

table_FC <- topTable(fit_contrasts, sort.by = "none", n = Inf) %>%
  select(c("P13vsP5",
           "P13vsP7"))

pv_native <- pathview(
  gene.data = table_FC,
  pathway.id = "hsa00532",
  species = "hsa",
  out.suffix = "P13vsP5andP7", 
  keys.align = "y",
  kegg.native = T,
  match.data = F,
  multi.state = T, 
  same.layer = F,
  low = list(gene = "dodgerblue")
)


pv_graphviz <- pathview(
  gene.data = table_FC,
  pathway.id = "hsa00532",
  species = "hsa",
  out.suffix = "P13vsP5andP7_graphviz", 
  keys.align = "y",
  kegg.native = F,
  match.data = F,
  multi.state = T, 
  same.layer = F,
  pdf.size = c(10, 10),
  low = list(gene = "dodgerblue")
)
