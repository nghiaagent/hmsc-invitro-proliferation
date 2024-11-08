# Install pacman, a package manager to install and load multiple packages
if (!requireNamespace("pacman")) {
  install.packages("pacman")
}

library(pacman)

# Install Bioconductor package manager
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

library(BiocManager)

# List Bioconductor packages to load/install
list_bioc_pkg <- c(
  "BiocHubsShiny",
  "AnnotationHub",
  "edgeR",
  "Homo.sapiens",
  "limma",
  "Glimma",
  "ComplexHeatmap",
  "tximeta",
  "SummarizedExperiment",
  "PCAtools",
  "sva",
  "ensembldb",
  "GSVA",
  "AnnotationDbi",
  "topGO",
  "pathview",
  "gage",
  "globaltest",
  "safe",
  "org.Hs.eg.db",
  "clusterProfiler",
  "ReactomePA",
  "reactome.db",
  "enrichplot",
  "GSEABase",
  "BiocStyle",
  "HDO.db",
  "DESeq2",
  "pcaExplorer",
  "IHW",
  "RUVSeq",
  "EnhancedVolcano"
)

# Install Bioconductor packages, if they are not yet installed
# Swap update = TRUE / FALSE depending on need to update pkg
BiocManager::install(
  pkgs = list_bioc_pkg,
  update = FALSE
)

# Load Bioconductor packages
invisible(lapply(list_bioc_pkg, function(x) {
  library(x, character.only = TRUE)
}))

# Load all CRAN packages
p_load(
  tibble,
  metan,
  tidyverse,
  magrittr,
  ggrepel,
  ggpubr,
  rstatix,
  ggalt,
  ggplotify,
  cowplot,
  gridExtra,
  statmod,
  volcano3D,
  msigdbr,
  viridis,
  palettetown,
  seriation,
  circlize,
  magick,
  cluster,
  rmdformats,
  shinybusy,
  WGCNA,
  data.table,
  DT,
  UpSetR,
  reactable,
  here,
  ashr,
  openxlsx,
  ggsignif
)

# Source other function or list definition scripts

source(here(
  "scripts",
  "0_define_GOIs.R"
))

# Clean up package list
rm(list_bioc_pkg)

# Limit number of cores used due to memory issues on laptops.
BiocParallel::register(SnowParam(workers = 4),
  default = TRUE
)
bpparam()

# Define contrasts for DESeq2
list_contrasts_deseq2 <- list(
  ## Coefs 1 - 6: Treatment at each timepoint
  Trt_P5_D3 = c("condition_ID", "P5D3Treated", "P5D3Untreated"),
  Trt_P5_D5 = c("condition_ID", "P5D5Treated", "P5D5Untreated"),
  Trt_P7_D3 = c("condition_ID", "P7D3Treated", "P7D3Untreated"),
  Trt_P7_D5 = c("condition_ID", "P7D5Treated", "P7D5Untreated"),
  Trt_P13_D3 = c("condition_ID", "P13D3Treated", "P13D3Untreated"),
  Trt_P13_D5 = c("condition_ID", "P13D5Treated", "P13D5Untreated"),

  ## Coefs 7 - 12: Day at each timepoint x treatment
  D5vsD3_UT_P5 = c("condition_ID", "P5D5Untreated", "P5D3Untreated"),
  D5vsD3_UT_P7 = c("condition_ID", "P7D5Untreated", "P7D3Untreated"),
  D5vsD3_UT_P13 = c("condition_ID", "P13D5Untreated", "P13D3Untreated"),
  D5vsD3_T_P5 = c("condition_ID", "P5D5Treated", "P5D3Treated"),
  D5vsD3_T_P7 = c("condition_ID", "P7D5Treated", "P7D3Treated"),
  D5vsD3_T_P13 = c("condition_ID", "P13D5Treated", "P13D3Treated"),

  ## Coefs 13 - 24: Passage at each day x treatment
  P7vsP5_UT_D3 = c("condition_ID", "P7D3Untreated", "P5D3Untreated"),
  P13vsP7_UT_D3 = c("condition_ID", "P13D3Untreated", "P7D3Untreated"),
  P13vsP5_UT_D3 = c("condition_ID", "P13D3Untreated", "P5D3Untreated"),
  P7vsP5_T_D3 = c("condition_ID", "P7D3Treated", "P5D3Treated"),
  P13vsP7_T_D3 = c("condition_ID", "P13D3Treated", "P7D3Treated"),
  P13vsP5_T_D3 = c("condition_ID", "P13D3Treated", "P5D3Treated"),
  P7vsP5_UT_D5 = c("condition_ID", "P7D5Untreated", "P5D5Untreated"),
  P13vsP7_UT_D5 = c("condition_ID", "P13D5Untreated", "P7D5Untreated"),
  P13vsP5_UT_D5 = c("condition_ID", "P13D5Untreated", "P5D5Untreated"),
  P7vsP5_T_D5 = c("condition_ID", "P7D5Treated", "P5D5Treated"),
  P13vsP7_T_D5 = c("condition_ID", "P13D5Treated", "P7D5Treated"),
  P13vsP5_T_D5 = c("condition_ID", "P13D5Treated", "P5D5Treated")
)

# Define functions to be used in the project

## Functions to format rowRanges and DESeqResults for required functions
format_deseq_results <- function(results) {
  # Type check
  if (class(results) != "DESeqResults") {
    stop("Provided results must be DESeqResults class")
  }

  # Format data, extract only relevant columns
  results_format <- results %>%
    as_tibble(rownames = "gene_id") %>%
    dplyr::select(c(gene_id, log2FoldChange, padj))

  # Return data
  return(results_format)
}

format_deseq_rowranges <- function(dds) {
  # Type check
  if (class(dds) != "DESeqDataSet") {
    stop("Provided dataset must be DESeqDataSet class")
  }

  # Format data, extract only relevant columns
  rowranges_format <- rowRanges(dds) %>%
    as_tibble() %>%
    dplyr::select(c(gene_id, gene_name, entrezid, description))

  # Return data
  return(rowranges_format)
}

## Function for extraction of top genes from a DESeqResults object
## sorted by LFC
## Can select top and bottom genes, or simply return all genes.
extract_topgenes <- function(
    results,
    dds,
    ntop = Inf,
    signif_only = TRUE) {
  # Format results and gene metadata
  results_format <- format_deseq_results(results)
  rowranges_format <- format_deseq_rowranges(dds)

  # Merge into one table, sort by decreasing fold change
  results_merge <- right_join(
    x = rowranges_format,
    y = results_format,
    by = join_by(gene_id == gene_id)
  ) %>%
    drop_na() %>%
    arrange(desc(log2FoldChange))

  # Remove non-significant genes
  if (signif_only == TRUE) {
    results_merge %<>% filter(padj < 0.05)
  }

  # Rename columns to be human readable
  results_merge %<>% dplyr::rename(
    "ENSEMBL ID" = gene_id,
    "Symbol" = gene_name,
    "Gene name" = description,
    "ENTREZ ID" = entrezid,
    "LogFC" = log2FoldChange,
    "adj. P-val" = padj
  )

  # Return split table (or not) depending on number of top genes provided
  # If not provided: Return entire table
  # If smaller number of top genes, return table containing top and bottom genes
  if (ntop == Inf) {
    results_merge %>%
      return()
  } else {
    results_merge %$%
      rbind(
        filter(., `LogFC` > 0) %>% slice_head(n = ntop),
        filter(., `LogFC` < 0) %>% slice_tail(n = ntop)
      ) %>%
      return()
  }
}

## Create function to join 2 DESeqResults from the list of DESeqResults
extract_joined_results <- function(
    results_1,
    results_2,
    results_lrt,
    name_1,
    name_2,
    dds) {
  # Format results and gene metadata
  rowranges_format <- format_deseq_rowranges(dds)

  results_1_format <- format_deseq_results(results_1)
  results_2_format <- format_deseq_results(results_2)

  results_lrt_format <- format_deseq_results(results_lrt) %>%
    dplyr::select(c(gene_id, padj)) %>%
    dplyr::rename(padj_lrt = padj)

  # Construct suffix for merged results
  suffix <- c(name_1, name_2) %>%
    map_chr(\(x) str_c("_", x))

  # Construct merged results table with inner join
  results_merge <- inner_join(
    x = results_1_format,
    y = results_2_format,
    by = join_by(gene_id == gene_id),
    suffix = suffix
  ) %>%
    right_join(
      x = results_lrt_format,
      y = .,
      by = join_by(gene_id == gene_id)
    ) %>%
    right_join(
      x = rowranges_format,
      y = .,
      by = join_by(gene_id == gene_id)
    ) %>%
    dplyr::rename(
      "ENSEMBL ID" = gene_id,
      "Symbol" = gene_name,
      "Gene name" = description,
      "ENTREZ ID" = entrezid,
    )

  # Return data
  return(results_merge)
}


## Define function to clip logFC and padj in DESeqResults to desired point
## For plotting with EnhancedVolcano
## Add relevant metadata
clip_results <- function(
    results,
    cutoff_logfc = 2.5,
    cutoff_padj = 1e-20,
    alpha = 0.05) {
  cutoff_logfc_neg <- cutoff_logfc * -1

  # Clip logFC to the threshold
  results$log2FoldChange %<>% case_when(
    . >= cutoff_logfc ~ cutoff_logfc,
    . <= cutoff_logfc_neg ~ cutoff_logfc_neg,
    .default = .
  )

  # Clip padj to the threshold
  results$padj %<>% case_when(
    . <= cutoff_padj ~ cutoff_padj,
    .default = .
  )

  # Add custom shapes to dots to identify clipped genes
  ## Normal dots: Shape 19
  ## Clipped (positive): Shape -9658
  ## Clipped (negative): Shape -9668
  ## Add names for legend
  results$volcano_shape <- case_when(
    results$log2FoldChange >= cutoff_logfc ~ -9658,
    results$log2FoldChange <= cutoff_logfc_neg ~ -9668,
    results$padj <= cutoff_padj ~ 17,
    .default = 19
  )

  names(results$volcano_shape) <- case_when(
    results$volcano_shape == -9658 ~ str_c("logFC >", cutoff_logfc),
    results$volcano_shape == -9668 ~ str_c("logFC < -", cutoff_logfc),
    results$volcano_shape == 17 ~ str_c("padj <", cutoff_padj),
    results$volcano_shape == 19 ~ "Unclipped"
  )

  # Make clipped genes larger
  results$volcano_size <- case_when(
    results$log2FoldChange >= cutoff_logfc ~ 3,
    results$log2FoldChange <= cutoff_logfc_neg ~ 3,
    results$padj <= cutoff_padj ~ 3,
    .default = 1
  )

  # Create new baseMean column for plotting with EnhancedVolcano
  results$baseMean_new <- 1 / (results$baseMean + 1)

  # Create new significance status level column
  results$colour <- case_when(
    results$padj <= alpha ~ "red2",
    .default = "grey30"
  )

  names(results$colour) <- case_when(
    results$colour == "red2" ~ str_c("padj <", alpha),
    results$colour == "grey30" ~ "ns"
  )

  # Return modified data
  return(results)
}
