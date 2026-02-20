### Don't source this file by itself;
### call in from another file after running env_prep.R
### Source
### Import salmon transcript quantification results
### to a single matrix, export that matrix.

# Construct sample table
table_samples <- here::here(
  "input",
  "annotation",
  "Cell_sample_table.csv"
) %>%
  read_csv() %>%
  mutate(
    ID = name,
    .keep = "unused"
  ) %>%
  mutate(
    ID = ID %>%
      str_replace(
        "hMSC_",
        "hMSC-"
      ) %>%
      str_replace(
        "(?<=0)_(?=[:digit:])",
        "-"
      )
  ) %>%
  mutate(
    cell_line = cell_line %>%
      factor(levels = c("hMSC-20176", "hMSC-21558")),
    Passage = Passage %>%
      factor(levels = c("P5", "P7", "P13")),
    Day = Day %>%
      factor(levels = c("D3", "D5")),
    Treatment = Treatment %>%
      factor(levels = c("Untreated", "Treated"))
  ) %>%
  mutate(
    run_date = run_date %>%
      str_replace_all(
        "_",
        ""
      )
  ) %>%
  mutate(
    timepoint_ID = str_c(Passage, Day, sep = "") %>%
      factor(
        levels = c(
          "P5D3",
          "P5D5",
          "P7D3",
          "P7D5",
          "P13D3",
          "P13D5"
        )
      )
  ) %>%
  mutate(
    condition_ID = str_c(Passage, Day, Treatment, sep = "") %>%
      factor(
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
      )
  ) %>%
  subset(included_in_dataset == TRUE) %>%
  arrange(ID, cell_line, Passage, Day, Treatment) %>%
  mutate(
    ID = ID %>%
      factor() %>%
      fct_inorder()
  ) %>%
  mutate(
    run_date = run_date %>%
      factor()
  )

# Construct transcriptome dataset for limma-voom
## Construct this for cDNA data from ENSEMBL
files_tx_quant <- table_samples$filename %>%
  str_c(
    "./input/cDNA/",
    .,
    "_quant/quant.sf"
  )

names_tx_quant <- files_tx_quant %>%
  str_replace(pattern = "IonXpress.*$", "") %>%
  str_replace(pattern = "./input/cDNA/", "") %>%
  str_replace(pattern = "_$", "")

list_files <- tibble(files = files_tx_quant, names = names_tx_quant) %>%
  cbind(dplyr::select(
    table_samples,
    c(
      "cell_line",
      "Passage",
      "Day",
      "Treatment",
      "run_date",
      "ID",
      "timepoint_ID",
      "condition_ID"
    )
  ))

all(file.exists(list_files$files)) # Make sure that all transcript files exist

# Import quantification results
# Summarise transcripts quantification to genes
# Create object for DGE with DESeq2
quant_deseq2 <- list_files %>%
  tximeta(type = "salmon") %>%
  summarizeToGene(assignRanges = "abundant") %>%
  DESeqDataSet(design = ~ condition_ID + run_date + cell_line)

# Filter lowly-expressed genes
# Keep only genes with higher than 10 counts in at least 3 samples
# Doesn't affect results; but speeds up computation
keep <- quant_deseq2 %$%
  filterByExpr(
    y = counts(.),
    group = colData(.)$condition_ID,
    min.count = 40,
    min.total.count = 60
  )

quant_deseq2_filter <- quant_deseq2[keep, ]

# Export quantification results
## Filtered
saveRDS(
  quant_deseq2_filter,
  file = here::here(
    "output",
    "data_expression",
    "pre_DGE",
    "quant_cDNA_deseq.RDS"
  )
)

## No filter
saveRDS(
  quant_deseq2,
  file = here::here(
    "output",
    "data_expression",
    "pre_DGE",
    "quant_cDNA_deseq_nofilter.RDS"
  )
)
