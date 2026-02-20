here::i_am("R/unused/1_pre_rowvars.R")

########################
# This script checks the variance of genes that are
# known to be good hMSC endogenous controls
# To see if we can apply RUV normalisation
########################

# Import packages
library(DESeq2)
library(here)
library(SummarizedExperiment)
library(tidyverse)


# Load dataset
# Remember that experimental design is already embedded in the dataset and model
# only need to extract comparisons
quant_deseq2 <- readRDS("output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")

# Define list of endogenous controls
list_ec <- c(
  "EEF1A1",
  "RPL13A",
  "RPLP0",
  "YWHAZ",
  "ACTB",
  "HPRT1",
  "GADD45A",
  "PUM1",
  "GAPDH",
  "TBP"
)

# Calculate variance of all genes
vars <- quant_deseq2 %>%
  counts() %>%
  rowVars()

vars_rlog <- quant_deseq2 %>%
  rlog() %>%
  assay() %>%
  rowVars()

## Combine into one table
gene_vars <- rowRanges(quant_deseq2) %>%
  as_tibble() %>%
  cbind(vars) %>%
  cbind(vars_rlog)

## Get variance ranks of endogenous controls
rank_by_counts <- list_ec %>%
  map(\(x) which(arrange(gene_vars, desc(vars))$symbol == x)) %>%
  set_names(list_ec)

rank_by_rlog <- list_ec %>%
  map(\(x) which(arrange(gene_vars, desc(vars_rlog))$symbol == x)) %>%
  set_names(list_ec)

# Save data
saveRDS(
  rank_by_counts,
  file = here::here(
    "output",
    "data_expression",
    "pre_DGE",
    "rank_by_counts.RDS"
  )
)

saveRDS(
  rank_by_rlog,
  file = here::here(
    "output",
    "data_expression",
    "pre_DGE",
    "rank_by_rlog.RDS"
  )
)
