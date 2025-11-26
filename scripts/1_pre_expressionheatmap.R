# Call DGE scripts first before calling this script
# Plots normalised log2CPM (quant_DGE_voom) as a heatmap.
# Use ?voom to check how it's made.
# limma authors recommend using logCPM from cpm() or voom()
# https://www.biostars.org/p/9511346/

source("./scripts/dge_cellpops_as_fixed.R")

siggenes <- topTable(fit_contrasts, number = Inf, sort.by = "F") %>%
    filter(adj.P.Val < 0.05)

quant_heatmap <- quant_DGE_voom[siggenes$GENEID, ]

E_heatmap <- t(scale(t(quant_heatmap$E)))

# Set color scheme and breaks

## For gene expression (z-scores)

col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

## For annotation data

col_cell_line <- palettetown::pokepal(6)[c(3, 4)]

col_treatment <- palettetown::pokepal(283)[c(11, 5)]

col_passage <- palettetown::pokepal(191)[c(8, 3, 9)]

col_day <- palettetown::pokepal(283)[c(3, 1)]

# Build annotation; include only necessary metadata

## Dataframe of annotation data

anno <- tibble(
    `Cell population` = quant_heatmap$targets$cell_line,
    `Passage` = quant_heatmap$targets$Passage,
    `Day` = quant_heatmap$targets$Day,
    `Treatment group` = quant_heatmap$targets$Treatment
)

## Colour mapping

anno_cols <- list(
    `Cell population` = c(
        "hMSC-20176" = col_cell_line[1],
        "hMSC-21558" = col_cell_line[2]
    ),
    `Passage` = c("P5" = col_passage[1], "P7" = col_passage[2], "P13" = col_passage[3]),
    `Day` = c("D3" = col_day[1], "D5" = col_day[2]),
    `Treatment` = c("Untreated" = col_treatment[1], "Treated" = col_treatment[2])
)

## ComplexHeatmap metadata annotation object

anno_object <- HeatmapAnnotation(
    df = anno,
    which = "col",
    col = anno_cols,
    annotation_height = 0.6,
    annotation_legend_param = list(
        `Cell population` = list(
            nrow = 2,
            title = "Cell population",
            title_position = "topleft",
            legend_direction = "vertical",
            title_gp = gpar(fontsize = 12, fontface = "bold"),
            labels_gp = gpar(fontsize = 12, fontface = "bold")
        ),
        `Passage` = list(
            nrow = 3,
            title = "Passage",
            title_position = "topleft",
            legend_direction = "vertical",
            title_gp = gpar(fontsize = 12, fontface = "bold"),
            labels_gp = gpar(fontsize = 12, fontface = "bold")
        ),
        `Day` = list(
            nrow = 2,
            title = "Day",
            title_position = "topleft",
            legend_direction = "vertical",
            title_gp = gpar(fontsize = 12, fontface = "bold"),
            labels_gp = gpar(fontsize = 12, fontface = "bold")
        ),
        `Treatment group` = list(
            nrow = 2,
            title = "Treatment",
            title_position = "topleft",
            legend_direction = "vertical",
            title_gp = gpar(fontsize = 12, fontface = "bold"),
            labels_gp = gpar(fontsize = 12, fontface = "bold")
        )
    )
)

# Select genes to be labelled

anno_genelabels <- rowAnnotation(
    Genes = anno_mark(
        at = seq(1, nrow(E_heatmap), 60),
        labels = rownames(E_heatmap)[seq(1, nrow(E_heatmap), 60)],
        labels_gp = gpar(fontsize = 8, fontface = "bold"),
        padding = 0.75
    ),
    width = unit(2.0, "cm") +

        max_text_width(rownames(E_heatmap)[seq(1, nrow(E_heatmap), 60)],
            gp = gpar(
                fontsize = 8, fontface = "bold"
            )
        )
)

# Apply row clustering

clusters <- pam(E_heatmap,
    k = 5
)

clusters$clustering <- paste0(
    "Cluster ",
    clusters$clustering
)

clusters$clustering <- factor(clusters$clustering,
    levels = c(
        "Cluster 1",
        "Cluster 2",
        "Cluster 3",
        "Cluster 4",
        "Cluster 5"
    )
)

# Create heatmap object

heatmap <- Heatmap(
    E_heatmap,
    name = "Gene\nZ-\nscore",
    col = colorRamp2(breaks, col),
    border = F,

    # Group rows based on clusters

    row_split = clusters$clustering,
    cluster_row_slices = FALSE,

    # parameters for the colour-bar that represents gradient of expression

    heatmap_legend_param = list(
        color_bar = "continuous",
        legend_direction = "vertical",
        legend_width = unit(8, "cm"),
        legend_height = unit(5.0, "cm"),
        title_position = "topcenter",
        title_gp = gpar(fontsize = 12, fontface = "bold"),
        labels_gp = gpar(fontsize = 12, fontface = "bold")
    ),

    # row (gene) parameters

    cluster_rows = T,
    show_row_dend = T,
    row_title = levels(clusters$clustering),
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_rot = 90,
    show_row_names = F,

    # column (sample) parameters

    cluster_column_slices = F,
    column_split = quant_heatmap$targets$Sample_ID,
    cluster_columns = T,
    show_column_dend = T,
    show_column_names = F,

    # specify top and bottom annotations

    top_annotation = anno_object
)

# Export heatmap


# Export heatmap of raw data

png(
    file = "./output/plots_heatmap/normalised_expression_genes_sample_ID.png",
    width = 12,
    height = 12,
    units = "in",
    res = 1200
)

export_heatmap <- plot(
    heatmap
)

dev.off()
