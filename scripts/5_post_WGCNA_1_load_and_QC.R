#### Run after DGE pipeline

# Transpose logCPM matrix for WGCNA
# Change from ENSEMBL to ENTREZ (to help with GSVA)

x <- quant_DGE_voom

rownames(x$genes) <- x$genes$GENEID

order <- order(fit$Amean, decreasing = TRUE)

x <- x[order,]

x <- x[complete.cases(dplyr::select(x$genes,!msigdb_h)),]

x <- x[!duplicated(x$genes$ENTREZID),]

rownames(x$E) <- x$genes$ENTREZID

rownames(x$genes) <- x$genes$ENTREZID

x <- x[order(x$genes$GENEID, decreasing = FALSE),]

x$genes <- x$genes %>% dplyr::select(!GENEID)

quant_DGE_WGCNA <- x

rm(x)

GCN_E <- t(quant_DGE_WGCNA$E) %>%
  as.data.frame()

# Check for low quality genes and samples
## Should be all OK
## If not, remove offending genes and samples from data

GCN_qual <- goodSamplesGenes(GCN_E)

if (!GCN_qual$allOK)
{
  # Print  genes and samples that are removed:
  if (sum(!GCN_qual$goodGenes) > 0)
    printFlush(paste("Removing genes:",
                     paste(names(datExpr0)[!GCN_qual$goodGenes],
                           collapse = ", ")))
  
  
  if (sum(!GCN_qual$goodSamples) > 0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!GCN_qual$goodSamples],
                           collapse = ", ")))
  
  
  # Remove the offending genes and samples from the data:
  GCN_E <- GCN_E[GCN_qual$goodSamples, GCN_qual$goodGenes]
}

# Perform sample clustering to detect outliers
## No outliers that are not already known (i.e. all related to batch effect)

meta_exclude <- c(
  'group',
  'norm.factors',
  'filename',
  'Sample_ID',
  'included_in_dataset',
  "lib.size",
  "ID",
  "timepoint_ID",
  "condition_ID"
)
png(
  "./output/plots_WGCNA/WGCNA_allsamples/Sample clustering.png",
  res = 150,
  width = 40,
  height = 25,
  units = "cm"
)

plotDendroAndColors(
  dendro = hclust(dist(GCN_E),
                  method = "average"),
  colors = labels2colors(dplyr::select(
    quant_DGE_WGCNA$targets,!all_of(meta_exclude)
  )),
  groupLabels = names(dplyr::select(
    quant_DGE_WGCNA$targets,!all_of(meta_exclude)
  )),
  main = "Hierarchical clustering of samples"
)

dev.off()

# Select soft thresholding power

## Calculate scale-free topology fit index

powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

GCN_sft <- pickSoftThreshold(
  GCN_E,
  powerVector = powers,
  verbose = 5
)

## Plot scale-free topology fit index and mean connectivity

png(
  "./output/plots_WGCNA/WGCNA_allsamples/Soft thresholding index selection.png",
  res = 150,
  width = 40,
  height = 20,
  units = "cm"
)

par(mfrow = c(1, 2))
cex1 = 0.9

### SFT index

plot(
  GCN_sft$fitIndices[, 1],-sign(GCN_sft$fitIndices[, 3]) * GCN_sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit,signed R^2",
  type = "n",
  main = paste("Scale independence")
)

text(
  GCN_sft$fitIndices[, 1],-sign(GCN_sft$fitIndices[, 3]) * GCN_sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)

abline(h = 0.90, col = "red")

### Mean connectivity

plot(
  GCN_sft$fitIndices[, 1],
  GCN_sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity")
)
text(
  GCN_sft$fitIndices[, 1],
  GCN_sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)

dev.off()

# I agree with the power selected in the object (power 7)

# Construct GCN object

GCN <- list(
  E = GCN_E,
  targets = quant_DGE_WGCNA$targets %>%
    mutate(
      cell_line = as.numeric(cell_line),
      Passage = as.numeric(Passage),
      Day = as.numeric(Day),
      Treatment = as.numeric(Treatment),
      run_date = as.numeric(factor(run_date)),
      timepoint_ID = as.numeric(timepoint_ID)
    ) %>%
    select(
      lib.size,
      cell_line,
      Passage,
      Day,
      Treatment,
      run_date,
      timepoint_ID
    ),
  genes = quant_DGE_WGCNA$genes,
  sft = GCN_sft
)

save(GCN,
     file = "./output/data_WGCNA/WGCNA_allsamples/GCN_input.Rdata")
