### Don't source this file by itself; call in from another file after running env_prep.R
### Source
### Import salmon transcript quantification results to a single matrix, export that matrix.

# Construct sample table

table_samples <-
  read_csv("./input/annotation/sample_table.csv") %>%
  mutate(cell_line = factor(
    cell_line,
    levels = c("WT HepG2-C3A cells", "ATP7B-KO HepG2-C3A cells")
  )) %>%
  mutate(genotype = factor(genotype, levels = c("WT", "ATP7B-KO"))) %>%
  mutate(treatment1 = factor(treatment1, levels = c("Untreated", "Cu"))) %>%
  mutate(treatment2 = factor(
    treatment2,
    levels = c("Untreated", "D-penicilamine", "trientine")
  )) %>%
  mutate(treatment_combined = factor(
    treatment_combined,
    levels = c("Untreated", "Cu", "Cu_D-penicilamine", "Cu_trientine")
  )) %>%
  mutate(condition_ID = factor(
    str_c(genotype, treatment_combined, sep = "_"),
    levels = c(
      "WT_Untreated",
      "WT_Cu",
      "WT_Cu_D-penicilamine",
      "WT_Cu_trientine",
      "ATP7B-KO_Untreated",
      "ATP7B-KO_Cu",
      "ATP7B-KO_Cu_D-penicilamine",
      "ATP7B-KO_Cu_trientine"
    )
  )) %>%
  arrange(condition_ID,
          genotype,
          treatment1,
          treatment2,
          treatment_combined,
          cell_line) %>%
  mutate(ID = factor(ID) %>% fct_inorder())

# Construct transcriptome dataset for limma-voom

## List files
## Make sure that all transcript files exist

list_files <- tibble(
  files = table_samples$ID %>%
    str_c("./input/cDNA/", ., ".genes.results"),
  names = table_samples$ID
)

all(file.exists(list_files$files))

# Import quantification results

quant_cDNA_gene <-
  tximeta(list_files, type = "rsem")

# Create object for DGE with edgeR/limma

quant_cDNA_DGE <-
  DGEList(assays(quant_cDNA_gene)[["counts"]])

## Add sample data

quant_cDNA_DGE$samples <-
  cbind(quant_cDNA_DGE$samples, table_samples)

## Add annotation data

geneid <- rownames(quant_cDNA_DGE)

genes <-
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = geneid,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL",
                "ENTREZID",
                "GENENAME",
                "SYMBOL",
                "GENETYPE")
  ) %>%
  dplyr::distinct(ENSEMBL, .keep_all = TRUE)

### Add inclusion of genes in MSigDB h

msigdb_h  <-
  read.gmt("./input/genesets/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.entrez.gmt") %>%
  mutate(term = as.character(term)) %>%
  mutate(ENTREZID = as.character(gene)) %>%
  mutate(gene = NULL) %>%
  nest_by(ENTREZID) %>%
  rowwise() %>%
  mutate(data = list(deframe(data))) %>%
  mutate(data = str_flatten_comma(data)) %>%
  mutate(msigdb_h = data) %>%
  mutate(data = NULL)

genes <- left_join(genes, msigdb_h)

quant_cDNA_DGE$genes <- genes

# Export quantification results

saveRDS(quant_cDNA_DGE, file = "./output/data_expression/pre_DGE/quant_cDNA_DGE.RDS")

saveRDS(quant_cDNA_gene, file = "./output/data_expression/pre_DGE/quant_cDNA_gene.RDS")
