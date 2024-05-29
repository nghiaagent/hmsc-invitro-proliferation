#### Run after DGE pipeline

# Transpose logCPM matrix for WGCNA

GCN_E <- t(quant_DGE_voom$E) %>%
  as.data.frame()

# Check for low quality genes and samples
## Should be all OK
## If not, remove offending genes and samples from data

GCN_qual <- goodSamplesGenes(quant_DGE_WGCNA)

if (!quant_DGE_qual$allOK)
{
  # Print  genes and samples that are removed:
  if (sum(!GCN_qual$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(datExpr0)[!GCN_qual$goodGenes], 
                           collapse = ", ")
                     )
               );
  
  if (sum(!GCN_qual$goodSamples)>0) 
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr0)[!GCN_qual$goodSamples],
                           collapse = ", ")
                     )
               );
  
  # Remove the offending genes and samples from the data:
  GCN_E <- GCN_E[GCN_qual$goodSamples, GCN_qual$goodGenes]
}

# Extract sample metadata

GCN_samples <- quant_DGE_voom$targets


# Perform sample clustering to detect outliers
## No outliers that are not already known (i.e. all related to batch effect)

meta_exclude <- c('group',
                  'norm.factors',
                  'filename',
                  'Sample_ID',
                  'included_in_dataset',
                  "lib.size",
                  "ID",
                  "timepoint_ID",
                  "condition_ID")
png(
  "./output/plots_WGCNA/Sample clustering.png",
  res = 150,
  width = 40,
  height = 25,
  units = "cm"
)

plotDendroAndColors(
  dendro = hclust(dist(GCN_E),
                  method = "average"),
  colors = labels2colors(dplyr::select(GCN_samples,
                                       !all_of(meta_exclude))),
  groupLabels = names(dplyr::select(GCN_samples,
                                    !all_of(meta_exclude))),
  main = "Hierarchical clustering of samples"
)

dev.off()
