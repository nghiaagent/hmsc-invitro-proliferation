# Load dataset
# Remember that experimental design is already embedded in the dataset and model
# only need to extract comparisons

quant_deseq2 <- readRDS("output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")

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

vars <- counts(quant_deseq2) %>%
  rowVars()

vars_rlog <- quant_deseq2 %>%
  rlog() %>%
  assay() %>%
  rowVars()

gene_vars <- rowRanges(quant_deseq2) %>%
  as_tibble() %>%
  cbind(vars) %>%
  cbind(vars_rlog)

rank_by_counts <- map(
  list_ec,
  \(x) which(arrange(gene_vars, desc(vars))$symbol == x)
) %>%
  set_names(list_ec)

rank_by_rlog <- map(
  list_ec,
  \(x) which(arrange(gene_vars, desc(vars_rlog))$symbol == x)
) %>%
  set_names(list_ec)
