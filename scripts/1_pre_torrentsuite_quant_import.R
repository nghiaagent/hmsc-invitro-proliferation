# Import transcript quantification in the form of barcode matrices output from Torrent Suite


# Construct sample table
# Populate barcode column based on filename

table_samples <-
  read_csv("./input/annotation/Cell_sample_table.csv") %>%
  mutate(barcode = str_extract(filename, "(?<=IonXpress_).*?(?=_)") %>%
           str_c("IonXpress_", ., sep = "")) %>%
  mutate(ID = name, .keep = "unused") %>%
  mutate(ID = str_replace(ID, "hMSC_", "hMSC-")) %>%
  mutate(ID = str_replace(ID, "(?<=0)_(?=[:digit:])", "-")) %>%
  mutate(cell_line = factor(cell_line, levels = c("hMSC-20176", "hMSC-21558"))) %>%
  mutate(Passage = factor(Passage, levels = c("P5", "P7", "P13"))) %>%
  mutate(Day = factor(Day, levels = c("D3", "D5"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Untreated", "Treated"))) %>%
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
  arrange(cell_line, Passage, Day, Treatment) %>%
  mutate(ID = fct_inorder(ID)) %>%
  subset(included_in_dataset == TRUE) %>%
  mutate(run_date = factor(run_date))

# Construct transcriptome dataset for limma-voom
## Construct this for data from Torrent Suite

names_tx_quant <- table_samples$run_date %>% unique() %>% sort()

barcodes <- table_samples %>%
  arrange(barcode) %>%
  group_by(run_date) %>%
  summarise(barcode = list(barcode), ID = list(as.character(ID))) %>%
  mutate(barcode = imap(.x = barcode, .f = \(x, idx) set_names(x, ID[[idx]])))

list_files <- tibble(name = names_tx_quant,
                     path = str_c("./input/torrentsuite/", names_tx_quant, ".tsv")) %>%
  left_join(barcodes, by = join_by(name == run_date))

# Test conversion of IDs from file to something useful

# Import quantification results

files <- pmap(.l = list_files, .f = \(name, path, barcode, ID) {
  read_tsv(path) %>%
    select(any_of(c("Gene", barcode)))
}) %>%
  set_names(list_files$name)

mat_torrentsuite <- files %>%
  reduce(left_join) %>%
  column_to_rownames("Gene") %>%
  select(any_of(levels(table_samples$ID))) %>%
  as.matrix()

quant_torrentsuite <- DESeqDataSetFromMatrix(
  countData = mat_torrentsuite,
  colData = table_samples %>%
    subset(included_in_dataset == TRUE),
  design = ~ condition_ID + run_date + cell_line
)

saveRDS(quant_torrentsuite, file = "./output/data_expression/pre_DGE/quant_torrentsuite_deseq.RDS")
