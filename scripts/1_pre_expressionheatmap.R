# Call DGE scripts first before calling this script
# Plots normalised log2CPM (quant_DGE_voom) as a heatmap. 
# Use ?voom to check how it's made.
# limma authors recommend using logCPM from cpm() or voom()
# https://www.biostars.org/p/9511346/

siggenes <- topTable(fit_contrasts, number = Inf, sort.by = "F") %>%
  filter(adj.P.Val < 0.05)

fit_heatmap <- fit_contrasts[siggenes$ENSEMBL, ]
fit_heatmap$EList <- fit_heatmap$EList[siggenes$ENSEMBL, ]

E_heatmap <- t(scale(t(fit_heatmap$EList$E)))

# Set color scheme and breaks

## For gene expression (z-scores)

col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

## For annotation data

col_treatment <- palettetown::pokepal(283)[c(11, 5, 7, 8)]

col_genotype <- palettetown::pokepal(283)[c(3, 1)]

# Build annotation; include only necessary metadata

## Dataframe of annotation data

anno <- tibble(
  `Genotype` = fit_heatmap$targets$genotype,
  `Treatment` = fit_heatmap$targets$treatment_combined
)

## Colour mapping

anno_cols <- list(
  `Genotype` = c(
    'WT' = col_genotype[1],
    'ATP7B-KO' = col_genotype[2]
  ),
  `Treatment` = c('Untreated' = col_treatment[1],
                  'Cu' = col_treatment[2],
                  'Cu_D-penicilamine' = col_treatment[3],
                  'Cu_trientine' = col_treatment[4]
                  )
)

## ComplexHeatmap metadata annotation object

anno_object <- HeatmapAnnotation(
  df = anno,
  which = 'col',
  col = anno_cols,
  annotation_height = 0.6,
  annotation_legend_param = list(
    `Genotype` = list(
      nrow = 2,
      title = 'Genotype',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')
    ),
    `Treatment` = list(
      nrow = 4,
      title = 'Treatment',
      title_position = 'topleft',
      legend_direction = 'vertical',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'bold')
    )
  )
)

# Select genes to be labelled

anno_genelabels <- rowAnnotation(
  Genes = anno_mark(
    at = seq(1, nrow(E_heatmap), 60),
    labels = rownames(E_heatmap)[seq(1, nrow(E_heatmap), 60)],
    labels_gp = gpar(fontsize = 8, fontface = 'bold'),
    padding = 0.75
  ),
  width = unit(2.0, 'cm') +
    
    max_text_width(rownames(E_heatmap)[seq(1, nrow(E_heatmap), 60)],
                   gp = gpar(
                     fontsize = 8,  fontface = 'bold'
                   ))
)

# Apply row clustering 

clusters <- pam(E_heatmap,
                k = 5)

clusters$clustering <- paste0('Cluster ', 
                              clusters$clustering)

clusters$clustering <- factor(clusters$clustering,
                              levels = c(
                                'Cluster 1',
                                'Cluster 2',
                                'Cluster 3',
                                'Cluster 4',
                                'Cluster 5'
                              ))

# Create heatmap object

## Unclustered

heatmap_unclustered <- Heatmap(
  E_heatmap,
  name = 'Gene\nZ-\nscore',
  col = colorRamp2(breaks, col),
  border = F,
  
  # Group rows based on clusters
  
  row_split = clusters$clustering,
  cluster_row_slices = FALSE,
  
  # parameters for the colour-bar that represents gradient of expression
  
  heatmap_legend_param = list(
    color_bar = 'continuous',
    legend_direction = 'vertical',
    legend_width = unit(8, 'cm'),
    legend_height = unit(5.0, 'cm'),
    title_position = 'topcenter',
    title_gp = gpar(fontsize = 12, fontface = 'bold'),
    labels_gp = gpar(fontsize = 12, fontface = 'bold')
  ),
  
  # row (gene) parameters
  
  cluster_rows = T,
  show_row_dend = T,      
  row_title = levels(clusters$clustering),
  row_title_side = 'left',
  row_title_gp = gpar(fontsize = 10,  fontface = 'bold'),
  row_title_rot = 90,
  show_row_names = F,
  
  # column (sample) parameters

  cluster_column_slices = F,
  cluster_columns = F,
  show_column_dend = T,
  show_column_names = F,
  
  # specify top and bottom annotations
  
  top_annotation = anno_object
)

# Export heatmap

png(
  file = "./output/plots_heatmap/normalised_expression_genes_unclustered.png",
  width = 12,
  height = 12,
  units = "in",
  res = 1200
)

export_heatmap <- plot(
  heatmap_unclustered
)

dev.off()

## Clustered

heatmap_clustered <- Heatmap(
  E_heatmap,
  name = 'Gene\nZ-\nscore',
  col = colorRamp2(breaks, col),
  border = F,
  
  # Group rows based on clusters
  
  row_split = clusters$clustering,
  cluster_row_slices = FALSE,
  
  # parameters for the colour-bar that represents gradient of expression
  
  heatmap_legend_param = list(
    color_bar = 'continuous',
    legend_direction = 'vertical',
    legend_width = unit(8, 'cm'),
    legend_height = unit(5.0, 'cm'),
    title_position = 'topcenter',
    title_gp = gpar(fontsize = 12, fontface = 'bold'),
    labels_gp = gpar(fontsize = 12, fontface = 'bold')
  ),
  
  # row (gene) parameters
  
  cluster_rows = T,
  show_row_dend = T,      
  row_title = levels(clusters$clustering),
  row_title_side = 'left',
  row_title_gp = gpar(fontsize = 10,  fontface = 'bold'),
  row_title_rot = 90,
  show_row_names = F,
  
  # column (sample) parameters
  
  cluster_column_slices = T,
  cluster_columns = T,
  show_column_dend = T,
  show_column_names = F,
  
  # specify top and bottom annotations
  
  top_annotation = anno_object
)

# Export heatmap

png(
  file = "./output/plots_heatmap/normalised_expression_genes_clustered.png",
  width = 12,
  height = 12,
  units = "in",
  res = 1200
)

export_heatmap <- plot(
  heatmap_clustered
)

dev.off()

