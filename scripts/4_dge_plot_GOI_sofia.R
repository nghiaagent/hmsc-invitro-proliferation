# Load data (batch-corrected counts)

source(here::here(
  "scripts",
  "0_meta_define_GOIs.R"
))

quant_deseq2_batchcor <- readRDS(file = here::here(
  "output",
  "data_expression",
  "post_DGE",
  "quant_deseq2_batchcor.RDS"
))

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
  filter(entrezid %in% geneids_goi)

# Define counts plotting function
## For treatments

plot_goi_treat <- function(gene_id, gene_name) {
  gene_counts <- plotCounts(
    quant_small,
    gene = gene_id,
    intgroup = c("Passage", "Treatment"),
    returnData = TRUE
  )

  plot <- ggplot(
    gene_counts,
    aes(
      x = Treatment,
      y = count,
      fill = Treatment
    )
  ) +
    geom_boxplot() +
    geom_jitter() +
    scale_y_log10() +
    theme_minimal_hgrid() +
    facet_wrap(~Passage) +
    ggtitle(label = gene_name) +
    scale_x_discrete(labels = c(
      "Treated" = "Hep",
      "Untreated" = "Ctrl"
    ))

  return(plot)
}

## For passages

plot_goi_passage <- function(gene_id, gene_name) {
  gene_counts <- plotCounts(
    quant_small,
    gene = gene_id,
    intgroup = c("Passage"),
    returnData = TRUE
  )

  plot <- ggplot(
    gene_counts,
    aes(
      x = Passage,
      y = count,
      fill = Passage
    )
  ) +
    geom_boxplot() +
    geom_jitter() +
    scale_y_log10() +
    theme_minimal_hgrid() +
    ggtitle(label = gene_name)

  return(plot)
}

## Apply to draw plots, export

plots_goi_treat <- map2(
  .x = genes_sel$gene_id,
  .y = genes_sel$gene_name,
  \(x, y) plot_goi_treat(x, y),
  .progress = TRUE
) %>%
  set_names(genes_sel$gene_name)

plots_goi_passage <- map2(
  .x = genes_sel$gene_id,
  .y = genes_sel$gene_name,
  \(x, y) plot_goi_passage(x, y),
  .progress = TRUE
) %>%
  set_names(genes_sel$gene_name)

## Export plots

purrr::walk2(
  .x = plots_goi_treat,
  .y = names(plots_goi_treat),
  \(x, y) {
    ggsave(
      filename = str_c(y, ".png", sep = ""),
      plot = x,
      path = here::here(
        ".",
        "output",
        "plots_boxplot_logCPM",
        "gois_treatment"
      ),
      scale = 0.7,
      width = 8,
      height = 8,
      units = "in",
      dpi = 144
    )
  },
  .progress = TRUE
)

purrr::walk2(
  .x = plots_goi_passage,
  .y = names(plots_goi_passage),
  \(x, y) {
    ggsave(
      filename = str_c(y, ".png", sep = ""),
      plot = x,
      path = here::here(
        ".",
        "output",
        "plots_boxplot_logCPM",
        "gois_passage"
      ),
      scale = 0.7,
      width = 8,
      height = 8,
      units = "in",
      dpi = 144
    )
  },
  .progress = TRUE
)
