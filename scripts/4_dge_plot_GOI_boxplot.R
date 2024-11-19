# Load data (batch-corrected counts)

quant_deseq2_batchcor <- readRDS(file = here::here(
  "output",
  "data_expression",
  "post_DGE",
  "quant_deseq2_batchcor.RDS"
))

results <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2.RDS"
  )
)

# Select data for conditions

conditions_interest <- c(
  "P5D3Untreated",
  "P7D3Untreated",
  "P13D3Untreated",
  "P5D3Treated",
  "P7D3Treated",
  "P13D3Treated"
)

quant_small <- quant_deseq2_batchcor %$%
  .[, colData(.)$condition_ID %in% conditions_interest]

# Extract rownames of GOIs

genes_sel <- rowRanges(quant_small) %>%
  as.data.frame() %>%
  filter(entrezid %in% geneids_goi_limited)

# Readjust DESeq results to include only GOIs

results_subset <- map(
  results,
  \(x) {
    x <- x[rownames(x) %in% genes_sel$gene_id, ]
    x$padj <- p.adjust(x$pvalue, method = "BH")

    return(x)
  }
)

# Define counts plotting functions

## For treatments
plot_goi_treat <- function(gene_id, gene_name) {
  # Get table containing normalised gene counts
  gene_counts <- plotCounts(
    quant_small,
    gene = gene_id,
    intgroup = c(
      "Passage",
      "Treatment"
    ),
    returnData = TRUE
  )

  # Get position to start drawing signif bars
  y_position <- log10(max(gene_counts$count) * 1.1)

  # Define position of bars and numbers manually

  annotation_signif <- tibble(
    Passage = factor(
      c("P5", "P7", "P13"),
      levels = c("P5", "P7", "P13")
    ),
    start = c("Untreated", "Untreated", "Untreated"),
    end = c("Treated", "Treated", "Treated"),
    y = c(y_position, y_position, y_position),
    label = c(
      results_subset[[1]][gene_id, ]$padj,
      results_subset[[3]][gene_id, ]$padj,
      results_subset[[5]][gene_id, ]$padj
    ) %>%
      stars_pval()
  )

  # Plot data
  plot <- ggplot(
    gene_counts,
    aes(x = Treatment, y = count)
  ) +
    geom_boxplot(aes(color = Treatment)) +
    geom_jitter(aes(color = Treatment)) +
    geom_signif(
      data = annotation_signif,
      aes(
        xmin = start,
        xmax = end,
        annotations = label,
        y_position = y
      ),
      tip_length = c(0.005, 0.005),
      manual = TRUE,
      textsize = 3
    ) +
    scale_x_discrete(labels = c(
      "Treated" = "Hep",
      "Untreated" = "Ctrl"
    )) +
    facet_wrap(~Passage) +
    scale_y_log10(expand = expansion(0, 0.05)) +
    theme_classic() +
    ggtitle(label = gene_name) +
    ylab("Normalised counts")

  # Return object
  return(plot)
}

## For passages

plot_goi_passage <- function(gene_id, gene_name) {
  # Get table containing normalised gene counts
  gene_counts <- plotCounts(
    quant_small,
    gene = gene_id,
    intgroup = c("Passage"),
    returnData = TRUE
  )

  # Get position to start drawing signif bars
  y_position <- log10(max(gene_counts$count))

  # Plot data
  plot <- ggplot(
    gene_counts,
    aes(x = Passage, y = count)
  ) +
    geom_boxplot(aes(color = Passage)) +
    geom_jitter(aes(color = Passage)) +
    geom_signif(
      comparisons = list(
        c("P5", "P7"),
        c("P7", "P13"),
        c("P5", "P13")
      ),
      annotation = c(
        results_subset[[13]][gene_id, ]$padj,
        results_subset[[14]][gene_id, ]$padj,
        results_subset[[15]][gene_id, ]$padj
      ) %>%
        stars_pval(),
      y_position = c(y_position, y_position, y_position),
      textsize = 3,
      step_increase = 0.1
    ) +
    scale_y_log10(expand = expansion(0, 0.05)) +
    theme_classic() +
    ggtitle(label = gene_name) +
    ylab("Normalised counts")

  # Return object
  return(plot)
}

## Apply to draw plots, export

plots_goi_treat <- map2(
  genes_sel$gene_id,
  genes_sel$gene_name,
  \(x, y) plot_goi_treat(x, y),
  .progress = TRUE
) %>%
  set_names(genes_sel$gene_name)

plots_goi_passage <- map2(
  genes_sel$gene_id,
  genes_sel$gene_name,
  \(x, y) plot_goi_passage(x, y),
  .progress = TRUE
) %>%
  set_names(genes_sel$gene_name)

## Export plots

purrr::iwalk(
  plots_goi_treat,
  \(x, idx) {
    ggsave(
      filename = str_c(idx, ".png", sep = ""),
      plot = x,
      path = here::here(
        ".",
        "output",
        "plots_boxplot_counts",
        "gois_treatment"
      ),
      scale = 0.6,
      width = 8,
      height = 6,
      units = "in",
      dpi = 144
    )
  },
  .progress = TRUE
)

purrr::iwalk(
  plots_goi_passage,
  \(x, idx) {
    ggsave(
      filename = str_c(idx, ".png", sep = ""),
      plot = x,
      path = here::here(
        ".",
        "output",
        "plots_boxplot_counts",
        "gois_passage"
      ),
      scale = 0.6,
      width = 8,
      height = 6,
      units = "in",
      dpi = 144
    )
  },
  .progress = TRUE
)

# Save data

saveRDS(
  plots_goi_treat,
  file = here::here(
    "output",
    "plots_boxplot_counts",
    "plots_goi_treat.RDS"
  )
)

saveRDS(
  plots_goi_passage,
  file = here::here(
    "output",
    "plots_boxplot_counts",
    "plots_goi_passage.RDS"
  )
)
