#### Load packages
### Install pacman, a package manager to install and load multiple packages at once

if (!requireNamespace("pacman")) {
  install.packages("pacman")
}
library(pacman)

### Some packages below use Bioconductor

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(BiocManager)

### List Bioconductor packages to be required

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
  "enrichplot"
)

### Install Bioconductor packages, if they are not yet installed and/or not up-to-date

BiocManager::install(
  pkgs = list_Bioc_Pkg,
  update = F
)

### Install packages specifically for EGSEA

p_load(
  HTMLUtils,
  hwriter,
  ggplot2,
  gplots,
  stringi,
  metap,
  devtools)

### Load Bioconductor packages with pacman

invisible(lapply(list_Bioc_Pkg, function(x)
  library(x, character.only = TRUE)))


### Load all CRAN packages with pacman

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
  cluster
)

#### Create folders for input and outputs, if not already present
# Code below means "if test for directory presence returns FALSE, create the directory"

### Check and create ./output/ folder

if (!dir.exists("./output")) {
  dir.create("./output")
}

if (!dir.exists("./output/plots_QC")) {
  dir.create("./output/plots_QC")
}

if (!dir.exists("./output/plots_PCA_prebatchcorrection")) {
  dir.create("./output/plots_PCA_prebatchcorrection")
}


if (!dir.exists("./output/plots_PCA_postbatchcorrection")) {
  dir.create("./output/plots_PCA_postbatchcorrection")
}

if (!dir.exists("./output/QC")) {
  dir.create("./output/QC")
}

if (!dir.exists("./output/GSEA")) {
  dir.create("./output/GSEA")
}

### Check and create ./input/  folder and subfolders

if (!dir.exists("./input")) {
  dir.create("./input")
}

# Clean up package list

rm(list_Bioc_Pkg)

# Limit number of cores used for multicore processing due to memory issues on laptops.

BiocParallel::register(SnowParam(workers = 8),
                       default = T)
bpparam()
