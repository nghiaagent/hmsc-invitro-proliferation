# Load data
fit_gsva <- readRDS(
    here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "GSVA_results.RDS"
    )
)

quant_gsva <- readRDS(
    here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "quant_GSVA.RDS"
    )
)

# Select relevant collections (h, c2, GOBP)
fit_gsva_subset <- fit_gsva[c("h", "c2_cp", "GOBP")]
quant_gsva_subset <- quant_gsva[c("h", "c2_cp", "GOBP")]

# Get significant gene sets
genesets_significant <- fit_gsva_subset %>%
    map(\(fit) {
        top <- fit %>%
            topTable(
                coef = NULL,
                number = Inf,
                sort.by = "none"
            ) %>%
            filter(adj.P.Val < 0.05) %>%
            rownames()

        # Return data
        return(top)
    })

quant_gsva_heatmap <- map2(
    quant_gsva_subset,
    genesets_significant,
    \(quant, genesets) quant[genesets, ]
)

es_heatmap <- quant_gsva_heatmap %>%
    map(\(x) {
        x %>%
            exprs() %>%
            t() %>%
            scale() %>%
            t()
    })

# Set color scheme and breaks
## For gene expression (rlog z-scores)
col <- inferno(n = 100)
breaks <- seq(-2, 2, length.out = 100)

# Create extra factors for splitting columns
split_cell_line <- quant_gsva_subset[[1]]@phenoData$cell_line

split_cell_line_passage <- str_c(
    quant_gsva_subset[[1]]@phenoData$cell_line,
    quant_gsva_subset[[1]]@phenoData$Passage,
    sep = "_"
) %>%
    factor(levels = c(
        "hMSC-20176_P5",
        "hMSC-20176_P7",
        "hMSC-20176_P13",
        "hMSC-21558_P5",
        "hMSC-21558_P7",
        "hMSC-21558_P13"
    ))

order <- quant_gsva_subset[[1]]@phenoData %>%
    pData() %>%
    arrange(cell_line, Passage, Day, Treatment) %>%
    rownames()

# Build annotation; include only necessary metadata
## Dataframe of annotation data
anno <- quant_gsva_subset[[1]]@phenoData %$%
    tibble(
        `Cell population` = .$cell_line,
        `Passage` = .$Passage,
        `Day` = .$Day,
        `Treatment` = .$Treatment %>%
            case_match(
                "Untreated" ~ "Control",
                "Treated" ~ "Heparin"
            )
    )

## Build ComplexHeatmap annotation object
anno_object <- HeatmapAnnotation(
    df = anno,
    which = "col",
    col = palette_heatmap,
    annotation_height = 0.6,
    annotation_legend_param = list(
        `Cell population` = list(
            nrow = 2,
            title = "Cell population",
            title_position = "topleft",
            legend_direction = "vertical",
            title_gp = gpar(fontsize = 8, fontface = "bold"),
            labels_gp = gpar(fontsize = 8, fontface = "bold")
        ),
        `Passage` = list(
            nrow = 3,
            title = "Passage",
            title_position = "topleft",
            legend_direction = "vertical",
            title_gp = gpar(fontsize = 8, fontface = "bold"),
            labels_gp = gpar(fontsize = 8, fontface = "bold")
        ),
        `Day` = list(
            nrow = 2,
            title = "Day",
            title_position = "topleft",
            legend_direction = "vertical",
            title_gp = gpar(fontsize = 8, fontface = "bold"),
            labels_gp = gpar(fontsize = 8, fontface = "bold")
        ),
        `Treatment` = list(
            nrow = 2,
            title = "Treatment",
            title_position = "topleft",
            legend_direction = "vertical",
            title_gp = gpar(fontsize = 8, fontface = "bold"),
            labels_gp = gpar(fontsize = 8, fontface = "bold")
        )
    )
)

# Create a heatmap to extract order for columns
heatmap_column_order <- Heatmap(
    es_heatmap[[3]],
    name = "ES\nZ-\nscore",
    col = colorRamp2(breaks, col),
    border = FALSE,

    # parameters for the colour-bar that represents gradient of expression
    heatmap_legend_param = list(
        color_bar = "continuous",
        legend_direction = "vertical",
        legend_width = unit(8, "cm"),
        legend_height = unit(5.0, "cm"),
        title_position = "topcenter",
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 8, fontface = "bold")
    ),

    # row (gene) parameters
    cluster_rows = TRUE,
    show_row_dend = FALSE,
    row_title_side = "left",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_title_rot = 90,
    show_row_names = FALSE,

    # column (sample) parameters
    column_split = split_cell_line,
    column_title = NULL,
    cluster_column_slices = FALSE,
    cluster_columns = TRUE,
    show_column_dend = TRUE,
    show_column_names = FALSE
) %>%
    plot() %>%
    column_order() %>%
    unlist()

# Create list of heatmap objects
list_heatmap <- es_heatmap %>%
    map(\(es) {
        heatmap <- Heatmap(
            es,
            name = "ES\nZ-\nscore",
            col = colorRamp2(breaks, col),
            border = FALSE,

            # parameters for the colour-bar that represents gradient of expression
            heatmap_legend_param = list(
                color_bar = "continuous",
                legend_direction = "vertical",
                legend_width = unit(8, "cm"),
                legend_height = unit(5.0, "cm"),
                title_position = "topcenter",
                title_gp = gpar(fontsize = 8, fontface = "bold"),
                labels_gp = gpar(fontsize = 8, fontface = "bold")
            ),

            # row (gene) parameters
            cluster_rows = TRUE,
            show_row_dend = FALSE,
            row_title_side = "left",
            row_title_gp = gpar(fontsize = 10, fontface = "bold"),
            row_title_rot = 90,
            show_row_names = FALSE,

            # column (sample) parameters
            column_split = split_cell_line_passage,
            column_order = order,
            column_title = NULL,
            cluster_column_slices = FALSE,
            cluster_columns = FALSE,
            show_column_dend = FALSE,
            show_column_names = FALSE,

            # specify top and bottom annotations
            top_annotation = anno_object
        )

        # Return data
        return(heatmap)
    })

# Export heatmaps
imap(
    list_heatmap,
    \(heatmap, name) {
        png(
            file = here::here(
                "output",
                "plots_heatmap",
                str_c("gsva_", name, "_split_cell_passage.png")
            ),
            width = 8,
            height = 4,
            units = "in",
            res = 600
        )

        plot(heatmap)

        dev.off()
    }
)
