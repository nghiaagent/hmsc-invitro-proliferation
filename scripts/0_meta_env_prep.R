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
  writexl,
  msigdbr,
  cowplot,
  viridis,
  palettetown,
  seriation,
  circlize,
  magick,
  cluster,
  rmdformats,
  future,
  shinybusy,
  WGCNA,
  data.table,
  DT,
  UpSetR,
  reactable,
  here,
  ashr
)

# Clean up package list

rm(list_bioc_pkg)

# Limit number of cores used due to memory issues on laptops.

BiocParallel::register(SnowParam(workers = 2),
  default = TRUE
)
bpparam()

# Define contrasts for DESeq2
## No interactions yet...

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
