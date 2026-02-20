here::i_am("R/4_dge_plot_03_quadrant.R")

########################
# Plot 4-way comparison plot to determine changes in gene expression
# between control and heparin-treated hMSCs over D3-D5
########################

# Import packages
library(conflicted)
library(cowplot)
library(DESeq2)
library(ggpmisc)
library(ggpp)
library(ggrepel)
library(here)
library(patchwork)
library(tidyverse)
library(tune)
library(writexl)

# Load data
results_lfcshrink <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2_lfcshrink.RDS"
  )
)

results_lrt <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_LRT.RDS"
  )
) %>%
  results(alpha = 0.05)

quant_deseq2_batchcor <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)
# Define coefficients to be included in plots
coef_results <- list(
  p5 = c(control = 7, heparin = 10),
  p7 = c(control = 8, heparin = 11),
  p13 = c(control = 9, heparin = 12)
)

# Create merged results tables
## Results for all genes
list_results <- coef_results %>%
  purrr::map(\(coefs) {
    extract_joined_results(
      results_1 = results_lfcshrink[[coefs[[1]]]],
      results_2 = results_lfcshrink[[coefs[[2]]]],
      results_lrt = results_lrt,
      name_1 = "control",
      name_2 = "heparin",
      quant_deseq2_batchcor
    ) %>%
      dplyr::mutate(
        outcome_control = dplyr::case_when(
          padj_control < 0.05 & log2FoldChange_control > 0 ~ 1,
          padj_control < 0.05 & log2FoldChange_control < 0 ~ -1,
          .default = 0
        ),
        outcome_heparin = dplyr::case_when(
          padj_heparin < 0.05 & log2FoldChange_heparin > 0 ~ 1,
          padj_heparin < 0.05 & log2FoldChange_heparin < 0 ~ -1,
          .default = 0
        )
      ) %>%
      dplyr::mutate(
        outcome_combined = dplyr::case_when(
          outcome_control != 0 & outcome_heparin == 0 ~ "control",
          outcome_control == 0 & outcome_heparin != 0 ~ "heparin",
          outcome_control != 0 & outcome_heparin != 0 ~ "both",
          outcome_control == 0 & outcome_heparin == 0 ~ "none"
        ) %>%
          factor(
            levels = c(
              "control",
              "heparin",
              "both",
              "none"
            ),
            labels = c(
              "Control only",
              "Heparin only",
              "Both",
              "ns"
            )
          )
      ) %>%
      dplyr::mutate(
        shape_combined = dplyr::case_when(
          outcome_combined == "ns" ~ 1,
          outcome_combined == "Both" ~ 2,
          outcome_combined == "Control only" & outcome_control == 1 ~ 3,
          outcome_combined == "Control only" & outcome_control == -1 ~ 4,
          outcome_combined == "Heparin only" & outcome_heparin == 1 ~ 5,
          outcome_combined == "Heparin only" & outcome_heparin == -1 ~ 6,
        ) %>%
          factor()
      )
  })

## Subset to all genes, genes within NF-kB, within TGF-beta
list_results_subset <- list_results %>%
  purrr::map(\(results) {
    results_subset <- list(
      all = results,
      nfkb = results %>% dplyr::filter(`ENTREZ ID` %in% geneids_nfkb),
      tgfb = results %>% dplyr::filter(`ENTREZ ID` %in% geneids_tgfb)
    )
  }) %>%
  unlist(recursive = FALSE)

# Create list of quadrant plots to
# show effect of Heparin over D3-D5 at each passage
list_plots <- purrr::map2(
  list_results,
  coef_results,
  \(results, coefs) {
    # Subset results to gene sets
    results_subset <- list(
      all = results,
      nfkb = results %>% dplyr::filter(`ENTREZ ID` %in% geneids_nfkb),
      tgfb = results %>% dplyr::filter(`ENTREZ ID` %in% geneids_tgfb)
    )
    plots <- results_subset %>%
      purrr::map(\(results) {
        results %>%
          dplyr::mutate(
            symbol_repel = dplyr::case_when(
              outcome_combined %in% c("Both", "ns") ~ NA,
              .default = Symbol
            )
          ) %>%
          ggplot2::ggplot(
            ggplot2::aes(
              x = log2FoldChange_control,
              y = log2FoldChange_heparin,
              fill = outcome_combined,
              label = symbol_repel
            )
          ) +
          ggplot2::geom_point(
            colour = "black",
            shape = 21,
            alpha = 0.7,
            stroke = 0.5,
            size = 2
          ) +
          ggrepel::geom_text_repel(
            ggplot2::aes(colour = outcome_combined),
            min.segment.length = 2,
            max.overlaps = 30,
            force = 4,
            force_pull = 0.3,
            size = 3,
            bg.color = "gray90",
            bg.r = 0.15
          ) +
          ggpp::geom_quadrant_lines(linetype = "dotted") +
          ggplot2::scale_fill_manual(values = palette_quadrant) +
          ggplot2::scale_colour_manual(values = palette_quadrant) +
          ggplot2::scale_x_continuous(limits = symmetric_limits(c(-2.6, 2.6))) +
          ggplot2::scale_y_continuous(limits = symmetric_limits(c(-2.6, 2.6))) +
          tune::coord_obs_pred() +
          ggplot2::theme_classic() +
          ggplot2::guides(shape = "none") +
          ggplot2::labs(
            x = str_c("log2 FC ", names(results_lfcshrink)[[coefs[[1]]]]),
            y = str_c("log2FC ", names(results_lfcshrink)[[coefs[[2]]]]),
            fill = "Outcome",
            colour = "Outcome"
          )
      })

    # Return data
    return(plots)
  }
) %>%
  unlist(recursive = FALSE)

# Build grid of plots and legend
grid_all <- list_plots %>%
  purrr::map(\(plot) plot + theme(legend.position = "none")) %>%
  patchwork::wrap_plots(ncol = 3)

grid_all <- (grid_all | cowplot::get_legend(list_plots[[1]])) +
  patchwork::plot_layout(widths = c(10, 1))

# Save data
## Gene list
writexl::write_xlsx(
  x = list_results_subset,
  path = here::here(
    "output",
    "data_expression",
    "genes_quadrant.xlsx"
  )
)

## Plot
ggsave(
  filename = here::here(
    "output",
    "plots_quadrant",
    "quadrant_all.png"
  ),
  height = 12,
  width = 12,
  scale = 0.9
)
