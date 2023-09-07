### Don't source this file by itself; call in from another file after running env_prep.R
### Source
### Import salmon transcript quantification results to a single matrix, export that matrix.

# Construct sample table

table_samples <-
  read_csv("./input/annotation/Cell_sample_table.csv") %>%
  mutate(ID = name,
         .keep = "unused") %>%
  dplyr::select(!c(name1)) %>%
  mutate(ID = str_replace(ID, "hMSC_", "hMSC-")) %>%
  mutate(ID = str_replace(ID, "(?<=0)_(?=[:digit:])", "-")) %>%
  mutate(cell_line = factor(cell_line,
                            levels = c("20176", "21558"))) %>%
  mutate(Passage = factor(Passage,
                          levels = c("P5", "P7", "P13"))) %>%
  mutate(Day = factor(Day,
                      levels = c("D3", "D5"))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Untreated", "Treated"))) %>%
  mutate(condition_ID = factor(
    str_c(Passage, Day, Treatment, sep = ""),
    levels = c(
      "P5D3Untreated",
      "P5D3Treated",
      "P5D5Untreated",
      "P5D5Treated",
      "P7D3Untreated",
      "P7D3Treated",
      "P7D5Untreated",
      "P7D5Treated",
      "P13D3Untreated",
      "P13D3Treated",
      "P13D5Untreated",
      "P13D5Treated",
    )
  )) %>%
  mutate(run_date = as_factor(run_date))

# Grab data from ENSEMBL109

ah <- AnnotationHub()
ensembl109 <- ah[["AH109606"]]

# Construct transcriptome dataset for limma-voom

## Construct this for cDNA data only

files_tx_quant <-
  list.files(
    path = "./input/cDNA",
    pattern = "*.sf",
    recursive = TRUE,
    full.names = T
  )

names_tx_quant <-
  str_replace(files_tx_quant, pattern = "IonXpress.*$", "") %>%
  str_replace(pattern = "./input/cDNA/", "") %>%
  str_replace(pattern = "_$", "")

list_files <- tibble(files = files_tx_quant,
                     names = names_tx_quant)

all(file.exists(list_files$files)) # Make sure that all transcript files exist

# Import quantification results

quant_cDNA_tx <-
  tximeta(list_files,
          type = "salmon")

# No issues

# Show dataset metadata
## Imported matrices
assayNames(quant_cDNA_tx)

## Transcripts imported
rowRanges(quant_cDNA_tx)

## Genomic information
seqinfo(quant_cDNA_tx)

## Database information
retrieveDb(quant_cDNA_tx) %>% class()

## The database originates from ENSEMBL

# Summarise transcripts quantification to genes

quant_cDNA_gene <- summarizeToGene(quant_cDNA_tx,
                                   countsFromAbundance = "lengthScaledTPM")

# Create object for DGE with edgeR
## Add sample data

quant_cDNA_DGE <- DGEList(assays(quant_cDNA_gene)[["counts"]])

quant_cDNA_DGE$samples$ID <- rownames(quant_cDNA_DGE$samples)

quant_cDNA_DGE$samples <-
  left_join(quant_cDNA_DGE$samples, table_samples)

## Add annotation data

geneid <- rownames(quant_cDNA_DGE)

genes <-
  AnnotationDbi::select(
    ensembl109,
    keys = geneid,
    keytype = "GENEID",
    columns = c(
      "ENTREZID",
      "GENENAME",
      "DESCRIPTION",
      "GENEBIOTYPE",
      "SEQNAME"
    )
  ) %>%
  dplyr::distinct(GENEID,
                  .keep_all = TRUE)

quant_cDNA_DGE$genes <- genes

# Export quantification results

saveRDS(quant_cDNA_DGE, file = "./output/quant_cDNA_DGE.RDS")

saveRDS(quant_cDNA_tx, file = "./output/quant_cDNA_tx.RDS")

saveRDS(quant_cDNA_gene, file = "./output/quant_cDNA_gene.RDS")

## Construct this for cDNA + ncRNA data from ENSEMBL

files_tx_quant <-
  list.files(
    path = "./input/cDNA_ncRNA_ENSEMBL/",
    pattern = "*.sf",
    recursive = TRUE,
    full.names = T
  )

names_tx_quant <-
  str_replace(files_tx_quant, pattern = "IonXpress.*$", "") %>%
  str_replace(pattern = "./input/cDNA_ncRNA_ENSEMBL/", "") %>%
  str_replace(pattern = "_$", "")

list_files <- tibble(files = files_tx_quant,
                     names = names_tx_quant)

all(file.exists(list_files$files)) # Make sure that all transcript files exist

# Import quantification results

quant_cDNA_ncRNA_ENSEMBL_tx <-
  tximeta(list_files,
          type = "salmon")

# No issues

# Show dataset metadata
## Imported matrices
assayNames(quant_cDNA_ncRNA_ENSEMBL_tx)

## Transcripts imported
rowRanges(quant_cDNA_ncRNA_ENSEMBL_tx)

## Genomic information
seqinfo(quant_cDNA_ncRNA_ENSEMBL_tx)

## Database information
retrieveDb(quant_cDNA_ncRNA_ENSEMBL_tx) %>% class()

## The database originates from ENSEMBL

# Summarise transcripts quantification to genes

quant_cDNA_ncRNA_ENSEMBL_gene <-
  summarizeToGene(quant_cDNA_ncRNA_ENSEMBL_tx,
                  countsFromAbundance = "lengthScaledTPM")

# Create object for DGE with edgeR
## Add sample data

quant_cDNA_ncRNA_ENSEMBL_DGE <-
  DGEList(assays(quant_cDNA_ncRNA_ENSEMBL_gene)[["counts"]])

quant_cDNA_ncRNA_ENSEMBL_DGE$samples$ID <-
  rownames(quant_cDNA_ncRNA_ENSEMBL_DGE$samples)

quant_cDNA_ncRNA_ENSEMBL_DGE$samples <-
  left_join(quant_cDNA_ncRNA_ENSEMBL_DGE$samples, table_samples)


## Add annotation data

geneid <- rownames(quant_cDNA_ncRNA_ENSEMBL_DGE)

genes <-
  AnnotationDbi::select(
    ensembl109,
    keys = geneid,
    keytype = "GENEID",
    columns = c(
      "ENTREZID",
      "GENENAME",
      "DESCRIPTION",
      "GENEBIOTYPE",
      "SEQNAME"
    )
  ) %>%
  dplyr::distinct(GENEID,
                  .keep_all = TRUE)

quant_cDNA_ncRNA_ENSEMBL_DGE$genes <- genes

# Export quantification results

saveRDS(quant_cDNA_ncRNA_ENSEMBL_DGE, file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE.RDS")

saveRDS(quant_cDNA_ncRNA_ENSEMBL_tx, file = "./output/quant_cDNA_ncRNA_ENSEMBL_tx.RDS")

saveRDS(quant_cDNA_ncRNA_ENSEMBL_gene, file = "./output/quant_cDNA_ncRNA_ENSEMBL_gene.RDS")
