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

### Load all packages with pacman

p_load(readr, dplyr, tidyr, tibble)

### List Bioconductor packages to be required

list_Bioc_Pkg <- c(
  "tximport",
  "edgeR",
  "Homo.sapiens",
  "limma",
  "ChAMP",
  "biomaRt",
  "goseq",
  "qusage",
  "maditr",
  "Glimma",
  "clusterProfiler",
  "EnsDb.Hsapiens.v75",
  "fgsea",
  "ComplexHeatmap"
)

### Install Bioconductor packages, if they are not yet installed and/or not up-to-date

BiocManager::install(
  pkgs = list_Bioc_Pkg,
  update = TRUE,
  ask = FALSE,
  force = FALSE
)

### Load Bioconductor packages with pacman

invisible(lapply(list_Bioc_Pkg, function(x)
  library(x, character.only = TRUE)))

#### Create folders for input and outputs, if not already present
# Code below means "if test for directory presence returns FALSE, create the directory"

### Check and create ./output/ folder

if (!dir.exists("./output")) {
  dir.create("./output")
}

### Check and create ./input/  folder

if (!dir.exists("./input")) {
  dir.create("./input")
}
