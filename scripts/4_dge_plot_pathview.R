# Prepare dataset. Run DGE first.
# Swap ENSEMBL ID with ENTREZ ID
# If multiple ENS genes share one ENTREZ ID, remove gene with lower average expression
# Remove genes with no Entrez ID

fit_contrasts <-
  readRDS(file = "./output/data_expression/post_DGE/fit_contrasts.RDS")

order <- order(fit_contrasts$Amean, decreasing = TRUE)

table_logFC <- topTable(
  fit_contrasts,
  coef = c(13, 14, 15),
  sort.by = "none",
  n = Inf
)

table_logFC <- table_logFC[order, ]

table_logFC <- table_logFC[complete.cases(table_logFC$ENTREZID), ]

table_logFC <- table_logFC[!duplicated(table_logFC$ENTREZID), ]

rownames(table_logFC) <- table_logFC$ENTREZID

table_logFC <- select(table_logFC,
                      c(P7vsP5_UT_D3,
                        P13vsP7_UT_D3,
                        P13vsP5_UT_D3))
# Draw pathway map

## Jump to output dir

setwd('./output/plots_kegg_pathways')

## Draw
### Set fun

plot_pv <- function(gene.data,
                    out.suffix) {
  pathview(
    gene.data = gene.data,
    pathway.id = c("hsa04064",
                   "hsa04668",
                   "hsa04010",
                   "hsa04310",
                   "hsa04350",
                   "hsa04110",
                   "hsa05205",
                   "hsa04512"
                   ),
    species = "hsa",
    kegg.dir = '../../input/kegg_pathways',
    kegg.native = TRUE,
    out.suffix = out.suffix,
    keys.align = "y",
    multi.state = TRUE,
    same.layer = TRUE,
    new.signature = FALSE,
    low = list(gene = "dodgerblue"),
    limit = list(gene = 2)
  )
}

### OPTIONAl - Draw 1 pair only. Code disabled for now
# 
# plot_pv(dplyr::select(table_logFC,
#                       1),
#         'P7vsP5')
# 
# plot_pv(dplyr::select(table_logFC,
#                       2),
#         'P13vsP7')
#
### 3 pairs

plot_pv(table_logFC,
        'P13vsP7vsP5')

## Return to wd

setwd('../../')
