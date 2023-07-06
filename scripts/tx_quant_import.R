### Don't source this file by itself; call in from another file after running env_prep.R
### Source 
### Import salmon transcript quantification results to a single matrix, export that matrix.

# Construct sample table

table_samples <- read_csv("./input/annotation/Cell_sample_table.csv") %>%
  mutate(ID = name,
         .keep = "unused") %>%
  dplyr::select(!c(name1)) %>%
  mutate(ID = str_replace(ID, "hMSC_", "hMSC-")) %>%
  mutate(ID = str_replace(ID, "(?<=0)_(?=[:digit:])", "-"))

# Grab data from ENSEMBL109

ensembl109 <- useMart(host="https://oct2022.archive.ensembl.org",
                      biomart="ENSEMBL_MART_ENSEMBL",
                      dataset="hsapiens_gene_ensembl") 

# Construct transcriptome dataset for limma-voom

## Construct this for cDNA data only

files_tx_quant <-
  list.files(path = "./input/cDNA",
             pattern = "*.sf",
             recursive = TRUE,
             full.names = T)

names_tx_quant <- str_replace(files_tx_quant, pattern = "IonXpress.*$", "") %>% 
  str_replace(pattern = "./input/cDNA/", "") %>% 
  str_replace(pattern = "_$", "")

list_files <- tibble(files = files_tx_quant,
                         names = names_tx_quant)

all(file.exists(list_files$files)) # Make sure that all transcript files exist

# Import quantification results

quant_cDNA_tx <-
  tximeta(
    list_files,
    type = "salmon"
  ) 

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
  
quant_cDNA_DGE$samples <- left_join(quant_cDNA_DGE$samples, table_samples)

## Add annotation data

geneid <- rownames(quant_cDNA_DGE)

genes <- select(ensembl109,
                keys=geneid,
                keytype="ensembl_gene_id", 
                columns=c("ensembl_gene_id", "external_gene_name", "description","chromosome_name","gene_biotype")) 

colnames(genes) <- c("ENSEMBL","SYMBOL","GENENAME","TXCHROM", "BIOTYPE")

quant_cDNA_DGE$genes <- genes

# Export quantification results

saveRDS(quant_cDNA_DGE, file = "./output/quant_cDNA_DGE.RDS")

saveRDS(quant_cDNA_tx, file = "./output/quant_cDNA_tx.RDS")

saveRDS(quant_cDNA_gene, file = "./output/quant_cDNA_gene.RDS")

## Construct this for cDNA + ncRNA data from ENSEMBL

files_tx_quant <-
  list.files(path = "./input/cDNA_ncRNA_ENSEMBL/",
             pattern = "*.sf",
             recursive = TRUE,
             full.names = T)

names_tx_quant <- str_replace(files_tx_quant, pattern = "IonXpress.*$", "") %>% 
  str_replace(pattern = "./input/cDNA_ncRNA_ENSEMBL/", "") %>% 
  str_replace(pattern = "_$", "")

list_files <- tibble(files = files_tx_quant,
                     names = names_tx_quant)

all(file.exists(list_files$files)) # Make sure that all transcript files exist

# Import quantification results

quant_cDNA_ncRNA_ENSEMBL_tx <-
  tximeta(
    list_files,
    type = "salmon"
  ) 

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

quant_cDNA_ncRNA_ENSEMBL_gene <- summarizeToGene(quant_cDNA_ncRNA_ENSEMBL_tx,
                                                 countsFromAbundance = "lengthScaledTPM")

# Create object for DGE with edgeR
## Add sample data

quant_cDNA_ncRNA_ENSEMBL_DGE <- DGEList(assays(quant_cDNA_ncRNA_ENSEMBL_gene)[["counts"]])

quant_cDNA_ncRNA_ENSEMBL_DGE$samples$ID <- rownames(quant_cDNA_ncRNA_ENSEMBL_DGE$samples)

quant_cDNA_ncRNA_ENSEMBL_DGE$samples <- left_join(quant_cDNA_ncRNA_ENSEMBL_DGE$samples, table_samples)


## Add annotation data

geneid <- rownames(quant_cDNA_ncRNA_ENSEMBL_DGE)

genes <- select(ensembl109,
                keys=geneid,
                keytype="ensembl_gene_id", 
                columns=c("ensembl_gene_id", "external_gene_name", "description","chromosome_name","gene_biotype")) 

colnames(genes) <- c("ENSEMBL","SYMBOL","GENENAME","TXCHROM", "BIOTYPE")

quant_cDNA_ncRNA_ENSEMBL_DGE$genes <- genes

# Export quantification results

saveRDS(quant_cDNA_ncRNA_ENSEMBL_DGE, file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE.RDS")

saveRDS(quant_cDNA_ncRNA_ENSEMBL_tx, file = "./output/quant_cDNA_ncRNA_ENSEMBL_tx.RDS")

saveRDS(quant_cDNA_ncRNA_ENSEMBL_gene, file = "./output/quant_cDNA_ncRNA_ENSEMBL_gene.RDS")

