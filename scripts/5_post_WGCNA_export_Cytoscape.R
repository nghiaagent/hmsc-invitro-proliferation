# Load data

## GCN with all samples

load(file.path(
  '.',
  'output',
  'data_WGCNA',
  'WGCNA_allsamples',
  'GCN_output.Rdata'
))

# Recalculate TOM

tom <- TOMsimilarityFromExpr(
  datExpr = GCN$E,
  power = GCN$sft$powerEstimate,
  TOMType = "unsigned",
  verbose = Inf
)

## Extract module 1 (Turquoise) nodes and edges

sel_turquoise <- labels2colors(GCN$net$colors) == "turquoise"

tom_turquoise <- tom[sel_turquoise, sel_turquoise]

dimnames(tom_turquoise) <- list(GCN$genes$GENENAME[sel_turquoise],
                                GCN$genes$GENENAME[sel_turquoise])

## Extract module 6 (Red) nodes and edges

sel_red <- labels2colors(GCN$net$colors) == "red"

tom_red <- tom[sel_red, sel_red]

dimnames(tom_red) <- list(GCN$genes$GENENAME[sel_red],
                          GCN$genes$GENENAME[sel_red])

## Export TOMs into edge and node lists for Cytoscape

### Turquoise net

exportNetworkToCytoscape(
  adjMat = tom_turquoise,
  edgeFile = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "cytoscape_input_edges_module_turquoise.txt"
  ),
  nodeFile = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "cytoscape_input_nodes_module_turquoise.txt"
  ),
  weighted = TRUE,
  threshold = 0.05,
  # Lower than normal threshold to allow filtering in Cytoscape
  nodeNames = GCN$genes$GENENAME[sel_turquoise],
  altNodeNames = GCN$genes$ENTREZID[sel_turquoise],
  nodeAttr = labels2colors(GCN$net$colors)[sel_turquoise]
)

### Red net

exportNetworkToCytoscape(
  adjMat = tom_red,
  edgeFile = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "cytoscape_input_edges_module_red.txt"
  ),
  nodeFile = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "cytoscape_input_nodes_module_red.txt"
  ),
  weighted = TRUE,
  threshold = 0.05,
  # Lower than normal threshold to allow filtering in Cytoscape
  nodeNames = GCN$genes$GENENAME[sel_red],
  altNodeNames = GCN$genes$ENTREZID[sel_red],
  nodeAttr = labels2colors(GCN$net$colors)[sel_red]
)
