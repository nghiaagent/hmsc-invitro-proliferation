### Don't source this file by itself; call in from another file after running env_prep.R

source("./scripts/dge_cellpops_as_fixed.R")

# Show list of comparisons

summary(decideTests(fit_contrasts))

# Extract list of genes for desired comparison - Phase C - A, D3

top_P13vsP5_UT_D3 <- topTable(fit_contrasts, coef = 5, n = Inf)

# Generate gene list for clusterProfiler. Nick ranked genes by t-statistic, I'll do the same.

## Extract lists of genes

### Genes with no ENTREZID are excluded from analyses that require ENTREZID

top_P13vsP5_UT_D3 <- dplyr::arrange(top_P13vsP5_UT_D3,
                                    desc(logFC))

list_gene_P13vsP5_UT_D3 <- top_P13vsP5_UT_D3$logFC

list_gene_P13vsP5_UT_D3_ENTREZ <- na.omit(top_P13vsP5_UT_D3)$logFC

names(list_gene_P13vsP5_UT_D3) <-
  as.character(top_P13vsP5_UT_D3$GENEID)

names(list_gene_P13vsP5_UT_D3_ENTREZ) <-
  as.character(na.omit(top_P13vsP5_UT_D3)$ENTREZID)

# Run GSEA with clusterProfiler

## Gene Ontology (GO) gene sets

GSEA_GOMF_P13vsP5_UT_D3 <- gseGO(
  list_gene_P13vsP5_UT_D3,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "MF",
  minGSSize = 25,
  maxGSSize = 500,
  pvalueCutoff = 0.05
) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()


GSEA_GOBP_P13vsP5_UT_D3 <- gseGO(
  list_gene_P13vsP5_UT_D3,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  minGSSize = 25,
  maxGSSize = 500,
  pvalueCutoff = 0.05
) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()


GSEA_GOCC_P13vsP5_UT_D3 <- gseGO(
  list_gene_P13vsP5_UT_D3,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "CC",
  minGSSize = 25,
  maxGSSize = 500,
  pvalueCutoff = 0.05
) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()

## KEGG gene sets

GSEA_KEGG_P13vsP5_UT_D3 <-
  gseKEGG(
    list_gene_P13vsP5_UT_D3_ENTREZ,
    organism = 'hsa',
    minGSSize = 25,
    pvalueCutoff = 0.05
  ) %>%
  setReadable('org.Hs.eg.db',
              'ENTREZID') %>%
  enrichplot::pairwise_termsim()

## ReactomePA gene sets

GSEA_ReactomePA_P13vsP5_UT_D3 <-
  gsePathway(list_gene_P13vsP5_UT_D3_ENTREZ,
             organism = 'human',
             minGSSize = 25) %>%
  setReadable('org.Hs.eg.db',
              'ENTREZID') %>%
  enrichplot::pairwise_termsim()

## WikiPathways gene sets

GSEA_WikiPathways_P13vsP5_UT_D3 <-
  gseWP(list_gene_P13vsP5_UT_D3_ENTREZ,
        organism = "Homo sapiens",
        minGSSize = 25) %>%
  setReadable('org.Hs.eg.db',
              'ENTREZID') %>%
  enrichplot::pairwise_termsim()

## MSigDB gene sets

### Obtain gene sets

#### set h: hallmark
#### set c2: curated
#### set c3: regulatory gene sets
#### set c5: ontology gene sets
#### set c8: cell type signature

msigdb_h <- msigdbr(species = "Homo sapiens",
                    category = "H") %>%
  dplyr::select(gs_name, ensembl_gene)

msigdb_c2 <- msigdbr(species = "Homo sapiens",
                     category = "C2") %>%
  dplyr::select(gs_name, ensembl_gene)

msigdb_c3 <- msigdbr(species = "Homo sapiens",
                     category = "C3") %>%
  dplyr::select(gs_name, ensembl_gene)

msigdb_c5 <- msigdbr(species = "Homo sapiens",
                     category = "C5") %>%
  dplyr::select(gs_name, ensembl_gene)

### Run GSEA

GSEA_MSigDB_h_P13vsP5_UT_D3 <-
  GSEA(list_gene_P13vsP5_UT_D3,
       TERM2GENE = msigdb_h,
       minGSSize = 25) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()

GSEA_MSigDB_c2_P13vsP5_UT_D3 <-
  GSEA(list_gene_P13vsP5_UT_D3,
       TERM2GENE = msigdb_c2,
       minGSSize = 25) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()

GSEA_MSigDB_c3_P13vsP5_UT_D3 <-
  GSEA(list_gene_P13vsP5_UT_D3,
       TERM2GENE = msigdb_c3,
       minGSSize = 25) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()

GSEA_MSigDB_c5_P13vsP5_UT_D3 <-
  GSEA(list_gene_P13vsP5_UT_D3,
       TERM2GENE = msigdb_c5,
       minGSSize = 25) %>%
  setReadable('org.Hs.eg.db',
              'ENSEMBL') %>%
  enrichplot::pairwise_termsim()

# Visualise results

dotplot(GSEA_GOBP_P13vsP5_UT_D3, showCategory = 30) + ggtitle("GO Biological Process \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_GOMF_P13vsP5_UT_D3, showCategory = 30) + ggtitle("GO Molecular Function \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_GOCC_P13vsP5_UT_D3, showCategory = 30) + ggtitle("GO Cellular Component \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_KEGG_P13vsP5_UT_D3, showCategory = 30) + ggtitle("KEGG \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_KEGG_P13vsP5_UT_D3, showCategory = 30) + ggtitle("KEGG \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_ReactomePA_P13vsP5_UT_D3, showCategory = 30) + ggtitle("ReactomePA \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_WikiPathways_P13vsP5_UT_D3, showCategory = 30) + ggtitle("WikiPathways \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_MSigDB_h_P13vsP5_UT_D3, showCategory = 30) + ggtitle("MSigDB Hallmarks \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_MSigDB_c2_P13vsP5_UT_D3, showCategory = 30) + ggtitle("MSigDB Curated Gene Sets \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_MSigDB_c3_P13vsP5_UT_D3, showCategory = 30) + ggtitle("MSigDB Regulatory Gene Sets \nPhase C - Phase A, Untreated, Day 3")

dotplot(GSEA_MSigDB_c5_P13vsP5_UT_D3, showCategory = 30) + ggtitle("MSigDB Ontology Gene Sets \nPhase C - Phase A, Untreated, Day 3")


# Pathways of interest - KEGG

hsa05022 <- pathview(gene.data  = list_gene_P13vsP5_UT_D3_ENTREZ,
                     pathway.id = "hsa05022",
                     species    = "hsa",
                     limit      = list(gene=max(abs(list_gene_P13vsP5_UT_D3_ENTREZ)), cpd=1))

hsa05010 <- pathview(gene.data  = list_gene_P13vsP5_UT_D3_ENTREZ,
                    pathway.id = "hsa05010",
                    species    = "hsa",
                    limit      = list(gene=max(abs(list_gene_P13vsP5_UT_D3_ENTREZ)), cpd=1))