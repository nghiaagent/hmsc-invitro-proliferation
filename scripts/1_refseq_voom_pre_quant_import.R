### Don't source this file by itself; call in from another file after running env_prep.R
### Source
### Import salmon transcript quantification results to a single matrix, export that matrix.

# Construct sample table

table_samples <-
  read_csv("./input/annotation/Cell_sample_table.csv") %>%
  mutate(ID = name, .keep = "unused") %>%
  mutate(ID = str_replace(ID, "hMSC_", "hMSC-")) %>%
  mutate(ID = str_replace(ID, "(?<=0)_(?=[:digit:])", "-")) %>%
  mutate(cell_line = factor(cell_line, levels = c("hMSC-20176", "hMSC-21558"))) %>%
  mutate(Passage = factor(Passage, levels = c("P5", "P7", "P13"))) %>%
  mutate(Day = factor(Day, levels = c("D3", "D5"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Untreated", "Treated"))) %>%
  mutate(run_date = str_replace_all(run_date, "_", "")) %>%
  mutate(timepoint_ID = factor(
    str_c(Passage, Day, sep = ""),
    levels = c("P5D3", "P5D5", "P7D3", "P7D5", "P13D3", "P13D5")
  )) %>%
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
      "P13D5Treated"
    )
  )) %>%
  subset(included_in_dataset == TRUE) %>%
  arrange(ID, cell_line, Passage, Day, Treatment) %>%
  mutate(ID = factor(ID)) %>%
  mutate(ID = fct_inorder(ID))

# Construct transcriptome dataset for limma-voom
## Construct this for cDNA + ncRNA data from ENSEMBL

files_tx_quant <- table_samples$filename

files_tx_quant <- str_c("./input/cDNA_refseq_filter/",
                        files_tx_quant,
                        "_quant/quant.sf")

names_tx_quant <-
  str_replace(files_tx_quant, pattern = "IonXpress.*$", "") %>%
  str_replace(pattern = "./input/cDNA_refseq_filter/", "") %>%
  str_replace(pattern = "_$", "")


list_files <- tibble(files = files_tx_quant, names = names_tx_quant)

all(file.exists(list_files$files)) # Make sure that all transcript files exist

# Prepare tx2gene list
# Maps a transcript to corresponding ENTREZID

keys <- read_tsv("./input/annotation/GCF_000001405.39_GRCh38.p13_rna_filter.fna.fai",
                 col_names = FALSE)$X1
tx2gene <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = keys,
  keytype = "REFSEQ",
  columns = c("REFSEQ", "ENTREZID")
)

# Import quantification results with tximport

quant <-
  tximport(
    list_files$files,
    type = "salmon",
    tx2gene = tx2gene,
    countsFromAbundance = "lengthScaledTPM"
  )

# Create object for DGE with edgeR
## Add sample data

quant_DGEList <- DGEList(quant[["counts"]])

quant_DGEList$samples <- cbind(quant_DGEList$samples, table_samples)

rownames(quant_DGEList$samples) <- quant_DGEList$samples$ID
colnames(quant_DGEList$counts) <- quant_DGEList$samples$ID

# Add annotation data

genes <- rownames(quant_DGEList) %>%
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = .,
    keytype = "ENTREZID",
    columns = c("ENTREZID", "GENENAME", "SYMBOL", "GENETYPE")
  ) %>%
  dplyr::distinct(ENTREZID, .keep_all = TRUE)

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

quant_DGEList$genes <- left_join(genes, msigdb_h)

# Export quantification results

saveRDS(quant_DGEList, file = "./output/data_expression/pre_DGE/quant_refseq_DGEList.RDS")