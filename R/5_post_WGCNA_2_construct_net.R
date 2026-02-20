here::i_am("R/5_post_WGCNA_2_construct_net.R")

########################
# Construct WGC network
########################

# Import packages
library(conflicted)
library(here)
library(tidyverse)
library(WGCNA)

# Set preferred implementations
## for cor()
conflicted::conflict_prefer("cor", "WGCNA")

# Load object
gcn <- readRDS(
  file = here::here(
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "GCN_input.Rdata"
  )
)

# Set desired power
gcn_power <- 8

# Construct GCN
gcn$net <- WGCNA::blockwiseModules(
  datExpr = gcn$E,
  power = gcn_power,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.35,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize = Inf,
  saveTOMs = TRUE,
  saveTOMFileBase = here::here(
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "blockwiseTOM"
  ),
  verbose = 5
)

# Plot dendrogram
png(
  here::here(
    "output",
    "plots_WGCNA",
    "WGCNA_allsamples",
    "GCN dendrogram.png"
  ),
  res = 150,
  width = 40,
  height = 20,
  units = "cm"
)

WGCNA::plotDendroAndColors(
  gcn$net$dendrograms[[1]],
  WGCNA::labels2colors(gcn$net$colors)[gcn$net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

dev.off()

# Correlate modules to samples metadata via heatmap plot
me_plot <- WGCNA::moduleEigengenes(
  gcn$E,
  WGCNA::labels2colors(gcn$net$colors)
) %>%
  .$eigengenes %>%
  WGCNA::orderMEs()

## Calculate pairwise Pearson cor of MEs against sample metadata
module2trait <- list(
  cor = WGCNA::cor(
    me_plot,
    gcn$targets,
    use = "p"
  ),
  pvalue = WGCNA::cor(
    me_plot,
    gcn$targets,
    use = "p"
  ) %>%
    WGCNA::corPvalueStudent(nrow(gcn$E))
)

## Display correlations and their p-values
png(
  here::here(
    "output",
    "plots_WGCNA",
    "WGCNA_allsamples",
    "Correlation of module PC1s to variables.png"
  ),
  res = 150,
  width = 30,
  height = 30,
  units = "cm"
)

text_matrix <- paste(
  signif(module2trait$cor, 2),
  "\n(",
  signif(module2trait$pvalue, 1),
  ")",
  sep = ""
)

dim(text_matrix) <- dim(module2trait$cor)
par(mar = c(6, 8.5, 3, 3), cex = 1.3)

## Display the correlation values within a heatmap plot
WGCNA::labeledHeatmap(
  Matrix = module2trait$cor,
  xLabels = names(gcn$targets),
  yLabels = names(me_plot),
  ySymbols = names(me_plot),
  colorLabels = FALSE,
  colors = WGCNA::blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-variable relationships")
)

dev.off()

# Save data
saveRDS(
  gcn,
  file = here::here(
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "GCN_output.RDS"
  )
)
