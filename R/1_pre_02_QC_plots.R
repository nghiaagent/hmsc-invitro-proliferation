here::i_am("R/1_pre_02_QC_plots.R")

########################
# This file creates and exports QC plots for the dataset
########################

# Import packages
library(data.table)
library(DEFormats)
library(DESeq2)
library(edgeR)
library(here)
library(SummarizedExperiment)
library(tidyverse)

# Load data
quant_cDNA_deseq <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "pre_DGE",
    "quant_cDNA_deseq.RDS"
  )
)
quant_cDNA_deseq_nofilter <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "pre_DGE",
    "quant_cDNA_deseq_nofilter.RDS"
  )
)

## Extract coldata
coldata <- colData(quant_cDNA_deseq) %>%
  as.data.frame()

# Get data for plotting libsize
## Construct a list of all possible normalisation methods
## Then turn list into a single dataframe
list_counts <- list(
  # DESeq2 median of ratios normalisation
  norm = quant_cDNA_deseq %>%
    DESeq2::estimateSizeFactors() %>%
    DESeq2::counts(normalized = TRUE),
  # Unnormalised data
  nonorm = quant_cDNA_deseq %>%
    DESeq2::counts(normalized = FALSE),
  # TMM normalisation
  norm_TMM = quant_cDNA_deseq %>%
    DEFormats::as.DGEList() %>%
    edgeR::normLibSizes() %>%
    edgeR::cpm(normalized.lib.sizes = FALSE),
  # Unfiltered DESeq2 median of ratios normalisation
  nofilter_norm = quant_cDNA_deseq_nofilter %>%
    DESeq2::estimateSizeFactors() %>%
    DESeq2::counts(normalized = TRUE),
  # Unfiltered unnormalised data
  nofilter_nonorm = quant_cDNA_deseq_nofilter %>%
    DESeq2::counts(normalized = FALSE)
) %>%
  # Turn list into dataframe for easier plotting
  imap(
    \(counts, name) {
      out <- counts %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "gene") %>%
        tidyr::pivot_longer(
          cols = !gene,
          names_to = "sample_id",
          values_to = "counts"
        ) %>%
        dplyr::mutate(logcounts_offset = log10(counts + 1)) %>%
        dplyr::mutate(data_type = name) %>%
        dplyr::left_join(
          coldata,
          by = dplyr::join_by("sample_id" == "names")
        )

      return(out)
    }
  ) %>%
  data.table::rbindlist()

# Calculate total library size
## Use unnormalised, unfiltered data only
df_libsize <- list_counts %>%
  dplyr::filter(data_type == "nofilter_nonorm") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(libsize = sum(counts)) %>%
  dplyr::left_join(
    coldata,
    by = dplyr::join_by("sample_id" == "names")
  ) %>%
  dplyr::arrange(cell_line) %>%
  dplyr::arrange(
    Passage,
    Day,
    Treatment
  ) %>%
  dplyr::mutate(
    sample_id = sample_id %>%
      factor() %>%
      forcats::fct_inorder()
  )

# Plot libsize
plot_libsize <- df_libsize %>%
  ggplot2::ggplot(ggplot2::aes(
    x = sample_id,
    y = libsize / 1000000,
    fill = timepoint_ID
  )) +
  ggplot2::geom_col(colour = "white") +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title.x = ggplot2::element_blank(),
    legend.title = ggplot2::element_text()
  ) +
  ggplot2::scale_fill_manual(values = palette_merge[1:6]) +
  ggplot2::geom_hline(
    yintercept = c(5, 10, 15),
    linetype = "dashed",
    color = "grey30"
  ) +
  ggplot2::labs(
    y = "Library size (million)",
    fill = "Timepoint"
  )

# Plot boxplots of normalised vs. unnormalised data
plots_normalisation <- list_counts %>%
  dplyr::filter(data_type %in% c("nonorm", "norm", "norm_TMM")) %>%
  dplyr::mutate(
    sample_id = sample_id %>%
      factor(levels = levels(df_libsize$sample_id))
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = sample_id,
      y = logcounts_offset,
      colour = timepoint_ID
    )
  ) +
  ggplot2::geom_boxplot() +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title.x = ggplot2::element_blank(),
    legend.title = ggplot2::element_text(),
    legend.position = "none"
  ) +
  ggplot2::scale_colour_manual(
    values = palette_merge[1:6]
  ) +
  ggplot2::labs(
    y = "log10(counts)"
  ) +
  ggplot2::facet_wrap(
    ~data_type,
    labeller = ggplot2::as_labeller(
      c(
        "nonorm" = "Unnormalised data",
        "norm" = "DESeq2 normalised data",
        "norm_TMM" = "TMM normalised data"
      )
    ),
    ncol = 1,
    scales = "free_y"
  )

# Plot density of genes before and after normalisation
plots_filtering <- list_counts %>%
  dplyr::filter(data_type %in% c("nofilter_norm", "norm")) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = logcounts_offset,
      colour = timepoint_ID
    )
  ) +
  ggplot2::geom_density(
    linewidth = 0.5,
    alpha = 0.6
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    legend.title = ggplot2::element_text(),
    legend.position = "none"
  ) +
  ggplot2::scale_colour_manual(
    values = palette_merge[1:6]
  ) +
  ggplot2::geom_vline(
    xintercept = log10(41),
    linetype = "dashed",
    colour = "grey30"
  ) +
  ggplot2::labs(
    x = "log10(counts)",
    y = "Density"
  ) +
  ggplot2::facet_wrap(
    ~data_type,
    ncol = 2,
    labeller = ggplot2::as_labeller(
      c(
        "nofilter_norm" = "Unfiltered data",
        "norm" = "FilterByExpr 40 counts"
      )
    ),
    scales = "free_y"
  )

# Save plots
## Library size
ggplot2::ggsave(
  here::here(
    "output",
    "plots_QC",
    "Library size.png"
  ),
  plot_libsize,
  width = 16,
  height = 9,
  scale = 0.7
)

## Distribution
ggplot2::ggsave(
  here::here(
    "output",
    "plots_QC",
    "Distribution.png"
  ),
  plots_normalisation,
  width = 8,
  height = 10,
  scale = 1
)

## Gene filtering
ggplot2::ggsave(
  here::here(
    "output",
    "plots_QC",
    "Gene filter.png"
  ),
  plots_filtering,
  width = 8,
  height = 6,
  scale = 0.7
)
