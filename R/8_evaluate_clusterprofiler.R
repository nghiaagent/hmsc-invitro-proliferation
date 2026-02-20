# Set variables
alpha <- 0.05
cutoff_logfc <- 10
cutoff_padj <- 1e-15

# Load data
results_deseq2 <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2.RDS"
  )
)

rowranges <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
) %>%
  rowRanges() %>%
  as.data.frame() %>%
  dplyr::select(
    gene_id,
    gene_name,
    entrezid
  )

contrasts_pilot <- c(
  "P7vsP5_UT_D3",
  "P13vsP7_UT_D3",
  "P13vsP5_UT_D3"
)

# Get desired contrasts, sort data
results_deseq2_pilot <- results_deseq2[contrasts_pilot] %>%
  map(\(.results) {
    .results <- .results %>%
      as.data.frame() %>%
      rownames_to_column(var = "gene_id") %>%
      as_tibble() %>%
      arrange(desc(stat)) %>%
      left_join(
        rowranges,
        by = join_by("gene_id" == "gene_id")
      )

    # Return data
    return(.results)
  })

## Create gene lists for ClusterProfiler
genelists_pilot <- results_deseq2_pilot %>%
  map(\(.results) {
    .results <- .results %>%
      dplyr::arrange(desc(stat))

    .gene_list <- .results$log2FoldChange %>%
      set_names(.results$gene_id)

    # Return data
    return(.gene_list)
  })

# Set up function for GSEA on GO gene sets
run_fgsea <- function(gene_list) {
  .gsea_gobp <- gseGO(
    gene_list,
    OrgDb = org.Hs.eg.db,
    keyType = "ENSEMBL",
    ont = "BP",
    minGSSize = 25,
    maxGSSize = 500,
    pvalueCutoff = 0.05
  ) %>%
    setReadable(
      "org.Hs.eg.db",
      "ENSEMBL"
    )

  # Return data
  return(.gsea_gobp)
}

# Run fgsea
results_fgsea <- genelists_pilot %>%
  map(\(.gene_list) run_fgsea(.gene_list))

# Create fgsea plots
plots_fgsea_dot <- results_fgsea %>%
  map(\(.fgsea) {
    if (nrow(.fgsea@result) == 0) {
      return(NA)
    } else {
      dotplot(.fgsea, showCategory = 15, label_format = 200)
    }
  })

plots_fgsea_cnet <- results_fgsea %>%
  map(\(.fgsea) {
    if (nrow(.fgsea@result) == 0) {
      return(NA)
    } else {
      cnetplot(.fgsea, showCategory = 15)
    }
  })

plots_fgsea_heat <- results_fgsea %>%
  map(\(.fgsea) {
    if (nrow(.fgsea@result) == 0) {
      return(NA)
    } else {
      heatplot(.fgsea, showCategory = 15)
    }
  })

plots_fgsea_worm <- results_fgsea %>%
  map(\(.fgsea) {
    if (nrow(.fgsea@result) == 0) {
      return(NA)
    } else {
      gseaplot2(
        .fgsea,
        geneSetID = c(
          "GO:0007608",
          "GO:0050906",
          "GO:0007606"
        )
      )
    }
  })

# Get volcano + MA plots of genes within sus sets
## Get geneids
gene_id_sus <- c(
  list_gmt[["GOBP"]][["GOBP_SENSORY_PERCEPTION_OF_SMELL"]]@geneIds,
  list_gmt[["GOBP"]][[
    "GOBP_DETECTION_OF_STIMULUS_INVOLVED_IN_SENSORY_PERCEPTION"
  ]]@geneIds,
  list_gmt[["GOBP"]][["GOBP_SENSORY_PERCEPTION_OF_CHEMICAL_STIMULUS"]]@geneIds
) %>%
  unique() %>%
  mapIds(
    org.Hs.eg.db,
    keys = .,
    column = "ENSEMBL",
    keytype = "ENTREZID"
  ) %>%
  .[!is.na(.)]

gene_id_moresus <- gene_id_sus %>%
  mapIds(
    org.Hs.eg.db,
    keys = .,
    column = "SYMBOL",
    keytype = "ENSEMBL"
  ) %>%
  str_which("^OR")

gene_id_sus <- gene_id_sus[gene_id_moresus]

## Clip results
results_deseq2_pilot_clipped <- results_deseq2[contrasts_pilot] %>%
  map(
    \(.results) {
      .results <- .results[complete.cases(.results), ]

      .results_clipped <- .results %>%
        clip_results(
          cutoff_logfc = cutoff_logfc,
          cutoff_padj = cutoff_padj,
          alpha = alpha
        )
    }
  )

plots_volcano_ma <- purrr::imap(
  results_deseq2_pilot_clipped,
  \(.results, .name) {
    .results <- .results[rownames(.results) %in% gene_id_sus, ]
    .labels <- rownames(.results) %>%
      mapIds(
        org.Hs.eg.db,
        keys = .,
        column = "SYMBOL",
        keytype = "ENSEMBL"
      )

    plots <- list(
      volcano = EnhancedVolcano(
        .results,
        x = "log2FoldChange",
        y = "padj",
        lab = .labels,
        xlim = c(cutoff_logfc * -1, cutoff_logfc),
        ylim = c(0, -log10(cutoff_padj)),
        ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
        axisLabSize = 8,
        title = .name,
        titleLabSize = 8,
        subtitle = NULL,
        caption = NULL,
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = .results$volcano_size,
        labSize = 2,
        boxedLabels = FALSE,
        shapeCustom = .results$volcano_shape,
        legendPosition = "none",
        drawConnectors = TRUE,
        widthConnectors = 0.1,
        colConnectors = "grey60",
        arrowheads = FALSE,
        gridlines.major = FALSE,
        gridlines.minor = FALSE
      ),
      ma = EnhancedVolcano(
        .results,
        x = "log2FoldChange",
        y = "baseMean_new",
        lab = .labels,
        selectLab = .labels[which(.results$padj <= alpha)],
        xlim = c(cutoff_logfc * -1, cutoff_logfc),
        ylim = c(0.5, 7),
        ylab = bquote(~ Log[10] ~ "mean of normalised counts"),
        axisLabSize = 8,
        title = .name,
        titleLabSize = 8,
        subtitle = NULL,
        caption = NULL,
        pCutoff = NA,
        FCcutoff = NA,
        pointSize = .results$volcano_size,
        labSize = 2,
        boxedLabels = FALSE,
        colCustom = .results$colour,
        shapeCustom = .results$volcano_shape,
        legendPosition = "none",
        drawConnectors = TRUE,
        widthConnectors = 0.1,
        colConnectors = "grey60",
        arrowheads = FALSE,
        gridlines.major = FALSE,
        gridlines.minor = FALSE
      )
    )

    # return data
    return(plots)
  },
  .progress = TRUE
)

# Construct plot grid
## First plot: Dot plot + GSEA worm demonstrating suspicious gene set
## Second plot: Dive into suspicious gene sets,
## demonstrating results may be invalid
plots_first <- wrap_plots(
  plots_fgsea_dot[[1]],
  plots_fgsea_dot[[2]],
  plots_fgsea_dot[[3]],
  plots_fgsea_worm[[1]][[1]],
  plots_fgsea_worm[[2]][[1]],
  plots_fgsea_worm[[3]][[1]],
  nrow = 3,
  byrow = FALSE
)

plots_second <- (plots_fgsea_cnet[[1]] |
  (wrap_plots(
    plots_volcano_ma %>%
      unlist(recursive = FALSE),
    nrow = 3
  ) +
    plot_layout(widths = c(2, 1)))) +
  plot_layout(widths = c(2, 4))


# Save plots
ggsave(
  here::here(
    "output",
    "plots_QC",
    "ClusterProfiler is bad 1.png"
  ),
  plots_first,
  width = 16,
  height = 9,
  scale = 1
)

ggsave(
  here::here(
    "output",
    "plots_QC",
    "ClusterProfiler is bad 2.png"
  ),
  plots_second,
  width = 16,
  height = 9,
  scale = 1
)
