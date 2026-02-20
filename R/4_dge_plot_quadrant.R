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
        "output", "data_expression", "post_DGE",
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
list_results <- coef_results %>%
    map(\(coefs) {
        extract_joined_results(
            results_1 = results_lfcshrink[[coefs[[1]]]],
            results_2 = results_lfcshrink[[coefs[[2]]]],
            results_lrt = results_lrt,
            name_1 = "control",
            name_2 = "heparin",
            quant_deseq2_batchcor
        ) %>%
            mutate(
                outcome_control = case_when(
                    padj_control < 0.05 & log2FoldChange_control > 0 ~ 1,
                    padj_control < 0.05 & log2FoldChange_control < 0 ~ -1,
                    .default = 0
                ),
                outcome_heparin = case_when(
                    padj_heparin < 0.05 & log2FoldChange_heparin > 0 ~ 1,
                    padj_heparin < 0.05 & log2FoldChange_heparin < 0 ~ -1,
                    .default = 0
                )
            ) %>%
            mutate(
                outcome_combined = case_when(
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
            mutate(
                shape_combined = case_when(
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

list_results_subset <- list_results %>%
    map(\(results) {
        results_subset <- list(
            all = results,
            nfkb = results %>% filter(`ENTREZ ID` %in% geneids_nfkb),
            tgfb = results %>% filter(`ENTREZ ID` %in% geneids_tgfb)
        )
    }) %>%
    unlist(recursive = FALSE)

# Create list of XY plots to show effect of Heparin over D3-D5 at each passage
list_plots <- map2(
    list_results,
    coef_results,
    \(results, coefs) {
        # Subset results to gene sets
        results_subset <- list(
            all = results,
            nfkb = results %>% filter(`ENTREZ ID` %in% geneids_nfkb),
            tgfb = results %>% filter(`ENTREZ ID` %in% geneids_tgfb)
        )
        plots <- results_subset %>%
            map(\(results) {
                results %>%
                    mutate(symbol_repel = case_when(
                        outcome_combined %in% c("Both", "ns") ~ NA,
                        .default = Symbol
                    )) %>%
                    ggplot(
                        aes(
                            x = log2FoldChange_control,
                            y = log2FoldChange_heparin,
                            fill = outcome_combined,
                            label = symbol_repel
                        )
                    ) +
                    geom_point(
                        colour = "black",
                        shape = 21,
                        alpha = 0.7,
                        stroke = 0.5,
                        size = 2
                    ) +
                    geom_text_repel(
                        aes(
                            colour = outcome_combined
                        ),
                        min.segment.length = 2,
                        max.overlaps = 30,
                        force = 4,
                        force_pull = 0.3,
                        size = 3,
                        bg.color = "gray90",
                        bg.r = 0.15
                    ) +
                    geom_quadrant_lines(linetype = "dotted") +
                    scale_fill_manual(values = palette_quadrant) +
                    scale_colour_manual(values = palette_quadrant) +
                    scale_x_continuous(limits = symmetric_limits(c(-2.6, 2.6))) +
                    scale_y_continuous(limits = symmetric_limits(c(-2.6, 2.6))) +
                    coord_obs_pred() +
                    theme_classic() +
                    guides(shape = "none") +
                    labs(
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

grid_all <- list_plots %>%
    map(\(plot) plot + theme(legend.position = "none")) %>%
    wrap_plots(ncol = 3)

grid_all <- (grid_all | get_legend(list_plots[[1]])) +
    plot_layout(widths = c(10, 1))

# Save data
## Gene list
write_xlsx(
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
