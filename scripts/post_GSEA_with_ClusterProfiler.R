### Don't source this file by itself; call in from another file after running env_prep.R

source("./scripts/dge_cellpops_as_fixed.R")

# Show list of comparisons

summary(decideTests(fit_contrasts))

# Extract list of genes for desired comparison - Phase C - A, D3

top_P13vsP5_UT_D3 <- topTable(fit_contrasts, coef = 5, n = Inf)

# Generate gene list for clusterProfiler. Nick ranked genes by t-statistic, I'll do the same.

top_P13vsP5_UT_D3 <- dplyr::arrange(top_P13vsP5_UT_D3,
                                    desc(logFC))

list_gene_P13vsP5_UT_D3 <- top_P13vsP5_UT_D3$logFC

names(list_gene_P13vsP5_UT_D3) <- as.character(top_P13vsP5_UT_D3$GENEID)

# Run GSEA with clusterProfiler

## Gene Ontology (GO) gene sets

GSEA_GOMF_P13vsP5_UT_D3 <- gseGO(list_gene_P13vsP5_UT_D3,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENSEMBL",
                               ont = "MF",
                               minGSSize = 100,
                               maxGSSize = 500,
                               pvalueCutoff = 0.05
                               ) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()


GSEA_GOBP_P13vsP5_UT_D3 <- gseGO(list_gene_P13vsP5_UT_D3,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENSEMBL",
                               ont = "BP",
                               minGSSize = 100,
                               maxGSSize = 500,
                               pvalueCutoff = 0.05
) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()


GSEA_GOCC_P13vsP5_UT_D3 <- gseGO(list_gene_P13vsP5_UT_D3,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENSEMBL",
                               ont = "CC",
                               minGSSize = 100,
                               maxGSSize = 500,
                               pvalueCutoff = 0.05
) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()

## KEGG gene sets

### Get Entrez ID for gene names

x <- bitr(names(list_gene_P13vsP5_UT_D3), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
x <- x[!duplicated(x$ENSEMBL),]
KEGG_mappings <- x[,2]
names(KEGG_mappings) <- x$ENSEMBL
rm(x)

GSEA_KEGG_P13vsP5_UT_D3 <- gseKEGG(list_gene_P13vsP5_UT_D3,
                                   organism = 'hsa',
                                   minGSSize = 100,
                                   pvalueCutoff = 0.05)

# Visualise results

## Results of GO

dotplot(GSEA_GOBP_P13vsP5_UT_D3, showCategory = 30) + ggtitle("GO Biological Process \nPhase C - Phase A, Untreated, D3")

dotplot(GSEA_GOMF_P13vsP5_UT_D3, showCategory = 30) + ggtitle("GO Molecular Function \nPhase C - Phase A, Untreated, D3")

dotplot(GSEA_GOCC_P13vsP5_UT_D3, showCategory = 30) + ggtitle("GO Cellular Component \nPhase C - Phase A, Untreated, D3")

emapplot(GSEA_GOBP_P13vsP5_UT_D3) + ggtitle("GO Biological Process \nPhase C - Phase A, Untreated, D3")

cnetplot(GSEA_GOBP_P13vsP5_UT_D3) + ggtitle("GO Biological Process \nPhase C - Phase A, Untreated, D3")

ridgeplot(GSEA_GOBP_P13vsP5_UT_D3) + ggtitle("GO Biological Process \nPhase C - Phase A, Untreated, D3")
