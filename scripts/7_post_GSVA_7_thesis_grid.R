# Load data
## Plots for passages
plots_poi_passage <- readRDS(
    here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "plots_poi_passage.RDS"
    )
) %>%
    unlist(recursive = FALSE)

# Create grid
## Get legend, reformat

plot_legend <- get_legend(
    plots_poi_passage[[1]] + guides(col = guide_legend(ncol = 3))
)

## Passages figure
plots_sel_poi_passage <- plots_poi_passage %>%
    map(\(plot) {
        plot +
            theme(
                legend.position = "none",
                plot.title = element_text(size = 8)
            )
    })


plots_sel_poi_passage <- wrap_plots(plots_sel_poi_passage) / plot_legend +
    plot_layout(heights = c(20, 1))

# Export plots

ggsave(
    filename = "pois_passage.png",
    plot = plots_sel_poi_passage,
    path = here::here(
        "output",
        "data_enrichment",
        "GSVA"
    ),
    scale = 0.8,
    width = 10,
    height = 11,
    units = "in",
    dpi = 144
)
