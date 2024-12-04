# Load data

## GCN with all samples
gcn <- readRDS(
  file = here::here(
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "GCN_output.RDS"
  )
)

# Recalculate TOM
tom <- TOMsimilarityFromExpr(
  datExpr = gcn$E,
  power = 8,
  TOMType = "unsigned",
  verbose = Inf
)

## Export TOMs into edge and node lists for Cytoscape
exportNetworkToCytoscape(
  adjMat = tom,
  edgeFile = here::here(
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "cytoscape_input_edges_module_all.txt"
  ),
  nodeFile = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "cytoscape_input_nodes_module_all.txt"
  ),
  weighted = TRUE,
  threshold = 0.1,
  # Lower than normal threshold to allow filtering in Cytoscape
  nodeNames = gcn$genes$symbol,
  altNodeNames = rownames(gcn$genes),
  nodeAttr = labels2colors(gcn$net$colors)
)
