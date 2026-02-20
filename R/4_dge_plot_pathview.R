here::i_am("R/4_dge_plot_pathview.R")

########################
# Load data
########################

# Import packages
library(DESeq2)
library(here)
library(tidyverse)

results_lfcshrink <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2_lfcshrink.RDS"
  )
)

quant_deseq2 <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

results_lrt <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_LRT.RDS"
  )
) %>%
  results(alpha = 0.05)

# Construct main table
results_betweenpassages <- extract_joined_results_trio(
  results_1 = results_lfcshrink[[13]],
  results_2 = results_lfcshrink[[14]],
  results_3 = results_lfcshrink[[15]],
  results_lrt = results_lrt,
  name_1 = "p5_p7",
  name_2 = "p7_p13",
  name_3 = "p5_p13",
  quant_deseq2
)

# Filter results table for use with Pathview
# Swap ENSEMBL ID with ENTREZ ID
# Remove gene with lower average expression (baseMean) if sharing ENSEMBL IDs
# Remove genes with no Entrez ID

order <- order(results_betweenpassages$baseMean, decreasing = TRUE)

results_betweenpassages <- results_betweenpassages %>%
  .[order, ] %>%
  .[map(.$`ENTREZ ID`, length) == 1, ] %>%
  .[!duplicated(.$`ENTREZ ID`), ] %>%
  as.data.frame()

results_pathview <- results_betweenpassages

rownames(results_pathview) <- results_pathview$`ENTREZ ID`

results_pathview <- results_pathview %>%
  dplyr::select(c(
    "log2FoldChange_p5_p7",
    "log2FoldChange_p7_p13",
    "log2FoldChange_p5_p13"
  ))

# Draw pathway map

## Jump to output dir
setwd("./output/plots_kegg_pathways")

## Draw
### Set fun

plot_pv <- function(gene_data, out_suffix) {
  pathview(
    gene.data = gene_data,
    pathway.id = c(
      "hsa04064",
      "hsa04668",
      "hsa04010",
      "hsa04310",
      "hsa04350",
      "hsa04110",
      "hsa05205",
      "hsa04512"
    ),
    species = "hsa",
    kegg.dir = "../../input/kegg_pathways",
    kegg.native = TRUE,
    out.suffix = out_suffix,
    keys.align = "y",
    multi.state = TRUE,
    same.layer = TRUE,
    new.signature = FALSE,
    low = list(gene = "dodgerblue"),
    limit = list(gene = 2)
  )
}

### 3 pairs

plot_pv(
  results_pathview,
  "P13vsP7vsP5"
)

## Return to wd

setwd("../../")
