# Install pacman, a package manager to install and load multiple packages at once

if (!requireNamespace("pacman")) {
  install.packages("pacman")
}
library(pacman)

# Install Bioconductor package manager

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocManager)

# List Bioconductor packages to load/install

list_Bioc_Pkg <- c(
  "BiocHubsShiny",
  "AnnotationHub",
  "edgeR",
  "Homo.sapiens",
  "limma",
  "biomaRt",
  "qusage",
  "maditr",
  "Glimma",
  "EnsDb.Hsapiens.v75",
  "ComplexHeatmap",
  "tximeta",
  "SummarizedExperiment",
  "PCAtools",
  "sva",
  "ensembldb",
  "PADOG",
  "GSVA",
  "AnnotationDbi",
  "topGO",
  "pathview",
  "gage",
  "globaltest",
  "safe",
  "org.Hs.eg.db",
  "variancePartition",
  "clusterProfiler",
  "ReactomePA",
  "enrichplot",
  "GSEABase",
  "BiocStyle"
)

# Install Bioconductor packages, if they are not yet installed
# Swap update = TRUE / FALSE depending on need to update pkg

BiocManager::install(
  pkgs = list_Bioc_Pkg,
  update = FALSE
)

# Load Bioconductor packages

invisible(lapply(list_Bioc_Pkg, function(x)
  library(x, character.only = TRUE)))

# Load all CRAN packages

p_load(
  tibble,
  tidyverse,
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
  shinybusy
)

# Clean up package list

rm(list_Bioc_Pkg)

# Limit number of cores used due to memory issues on laptops.

BiocParallel::register(SnowParam(workers = 8),
                       default = T)
bpparam()

# TO ADD: Create folders needed for outputting files
