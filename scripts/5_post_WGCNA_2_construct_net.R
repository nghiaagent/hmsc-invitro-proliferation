# Load object

load(file = "./output/data_WGCNA/WGCNA_allsamples/GCN_input.Rdata")

# Construct GCN

GCN$net <- blockwiseModules(
  datExpr = GCN$E,
  power = GCN$sft$powerEstimate,
  TOMType = "unsigned",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.35,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  maxBlockSize = Inf,
  saveTOMs = TRUE,
  saveTOMFileBase = file.path(
    '.',
    'output',
    'data_WGCNA',
    'WGCNA_allsamples',
    'blockwiseTOM'
  ),
  verbose = 5
)

# Plot dendrogram

png(
  "./output/plots_WGCNA/WGCNA_allsamples/GCN dendrogram.png",
  res = 150,
  width = 40,
  height = 20,
  units = "cm"
)

plotDendroAndColors(
  GCN$net$dendrograms[[1]],
  labels2colors(GCN$net$colors)[GCN$net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
)

dev.off()

# Correlate modules to samples metadata

MEs_plot <- moduleEigengenes(GCN$E, labels2colors(GCN$net$colors))$eigengenes %>%
  orderMEs()

## Calculate pairwise Pearson cor of MEs against sample metadata

module2trait <- list(
  cor = cor(MEs_plot, GCN$targets, use = "p"),
  pvalue = cor(MEs_plot, GCN$targets, use = "p") %>%
    corPvalueStudent(nrow(GCN$E))
)

# Will display correlations and their p-values

png(
  "./output/plots_WGCNA/WGCNA_allsamples/Correlation of module PC1s to variables.png",
  res = 150,
  width = 30,
  height = 30,
  units = "cm"
)

textMatrix =  paste(signif(module2trait$cor, 2),
                    "\n(",
                    signif(module2trait$pvalue, 1),
                    ")",
                    sep = "")

dim(textMatrix) = dim(module2trait$cor)
par(mar = c(6, 8.5, 3, 3), cex = 1.3)

# Display the correlation values within a heatmap plot

labeledHeatmap(
  Matrix = module2trait$cor,
  xLabels = names(GCN$targets),
  yLabels = names(MEs_plot),
  ySymbols = names(MEs_plot),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-variable relationships")
)

dev.off()

save(GCN, file = "./output/data_WGCNA/WGCNA_allsamples/GCN_output.RData")