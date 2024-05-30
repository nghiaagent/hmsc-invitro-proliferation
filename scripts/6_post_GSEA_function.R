# Run relevant scripts beforehand
# dge_cellpops_as_fixed
# dge_extract_ordered_genelist

# List of gene sets to be used for GSEA
# GO (CC, BP, MF)
# KEGG
# ReactomePA
# MSigDB h, c2, c3, c5

# Obtain gene sets for MSigDB

## set h: hallmark
## set c2: curated
## set c3: regulatory gene sets
## set c5: ontology gene sets

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

# Define function to run GSEA on (ordered) list of genes extracted from topTable
# Takes gene list as input
# Outputs folder in ./output/data_enrichment/GSEA with all enrichment RDS objects
# And gene lists ready for Cytoscape


run_GSEA <- function(genes_list,
                     name_output) {
  # Make prerequisite folders under ./output
  
  name_output <- as.character(name_output)
  
  path_output <- file.path(getwd(),
                           'output',
                           'data_enrichment',
                           'GSEA',
                           name_output)
  
  message(str_c("Output GSEA results to",
                path_output,
                sep = " "))
  
  if (!dir.exists(path_output)) {
    dir.create(path_output)
  }
  
  # Run GSEA on GO gene sets
  
  message(str_c("Running GOMF enrichment for",
                name_output,
                sep = " "))
  
  GSEA_GOMF <- gseGO(
    genes_list$GENEID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "MF",
    minGSSize = 25,
    maxGSSize = 500,
    pvalueCutoff = 0.05
  ) %>%
    setReadable('org.Hs.eg.db',
                'ENSEMBL')
  
  message(str_c("Running GOBP enrichment for",
                name_output,
                sep = " "))
  
  GSEA_GOBP <- gseGO(
    genes_list$GENEID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    minGSSize = 25,
    maxGSSize = 500,
    pvalueCutoff = 0.05
  ) %>%
    setReadable('org.Hs.eg.db',
                'ENSEMBL')
  
  message(str_c("Running GOCC enrichment for",
                name_output,
                sep = " "))
  
  GSEA_GOCC <- gseGO(
    genes_list$GENEID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "CC",
    minGSSize = 25,
    maxGSSize = 500,
    pvalueCutoff = 0.05
  ) %>%
    setReadable('org.Hs.eg.db',
                'ENSEMBL')
  
  # Run GSEA on KEGG
  
  message(str_c("Running KEGG enrichment for",
                name_output,
                sep = " "))
  
  GSEA_KEGG <- gseKEGG(
    genes_list$ENTREZID,
    organism = 'hsa',
    minGSSize = 25,
    pvalueCutoff = 0.05
  ) %>%
    setReadable('org.Hs.eg.db',
                'ENTREZID')
  
  # Run GSEA on ReactomePA
  
  message(str_c("Running Reactome enrichment for",
                name_output,
                sep = " "))
  
  GSEA_Reactome <-
    gsePathway(genes_list$ENTREZID,
               organism = 'human',
               minGSSize = 25) %>%
    setReadable('org.Hs.eg.db',
                'ENTREZID')
  
  # Run GSEA on MSigDB
  
  message(str_c("Running MSigDB h enrichment for",
                name_output,
                sep = " "))
  
  GSEA_MSigDB_h <-
    GSEA(genes_list$GENEID,
         TERM2GENE = msigdb_h,
         minGSSize = 25) %>%
    setReadable('org.Hs.eg.db',
                'ENSEMBL')
  
  message(str_c("Running MSigDB c2 enrichment for",
                name_output,
                sep = " "))
  
  GSEA_MSigDB_c2 <-
    GSEA(genes_list$GENEID,
         TERM2GENE = msigdb_c2,
         minGSSize = 25) %>%
    setReadable('org.Hs.eg.db',
                'ENSEMBL')
  
  message(str_c("Running MSigDB c3 enrichment for",
                name_output,
                sep = " "))
  
  GSEA_MSigDB_c3 <-
    GSEA(genes_list$GENEID,
         TERM2GENE = msigdb_c3,
         minGSSize = 25) %>%
    setReadable('org.Hs.eg.db',
                'ENSEMBL')
  
  message(str_c("Running MSigDB c5 enrichment for",
                name_output,
                sep = " "))
  
  GSEA_MSigDB_c5 <-
    GSEA(genes_list$GENEID,
         TERM2GENE = msigdb_c5,
         minGSSize = 25) %>%
    setReadable('org.Hs.eg.db',
                'ENSEMBL')
  
  # Export RDS
  
  message(str_c("Saving RDS to",
                name_output,
                sep = " "))
  
  saveRDS(
    list(
      "GOBP" = GSEA_GOBP,
      "GOCC" = GSEA_GOCC,
      "GOMF" = GSEA_GOMF,
      "KEGG" = GSEA_KEGG,
      "MSigDB_h" = GSEA_MSigDB_h,
      "MSigDB_c2" = GSEA_MSigDB_c2,
      "MSigDB_c3" = GSEA_MSigDB_c3,
      "MSigDB_c5" = GSEA_MSigDB_c5,
      "Reactome" = GSEA_Reactome,
      "path" = as.character(path_output)
    ),
    file = file.path(path_output,
                     "GSEA_results.RDS")
  )
}