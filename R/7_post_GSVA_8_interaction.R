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

# Set conditions of interest
conditions_interest <- c(
    ## Coefs 7 - 12: Day at each timepoint x treatment
    "D5vsD3_UT_P5" = 7,
    "D5vsD3_UT_P7" = 8,
    "D5vsD3_UT_P13" = 9,
    "D5vsD3_T_P5" = 10,
    "D5vsD3_T_P7" = 11,
    "D5vsD3_T_P13" = 12
)
names_conditions_interest <- list(
    p5 = c("D5vsD3_UT_P5", "D5vsD3_T_P5"),
    p7 = c("D5vsD3_UT_P7", "D5vsD3_T_P7"),
    p13 = c("D5vsD3_UT_P13", "D5vsD3_T_P13")
)

# Create nested list containing TopTables
toptables_gsva <- map(
    conditions_interest,
    \(coef) {
        topTable(
            fit_gsva[[1]],
            coef = coef,
            number = Inf,
            sort.by = "none"
        ) %>%
            rownames_to_column(var = "set")
    }
)

list_results <- map2(
    toptables_gsva[1:3],
    toptables_gsva[4:6],
    \(toptable_control, toptable_heparin) {
        full_join(
            toptable_control,
            toptable_heparin,
            by = join_by(set == set),
            suffix = c("_control", "_heparin")
        ) %>%
            mutate(
                outcome_control = case_when(
                    adj.P.Val_control < 0.05 & logFC_control > 0 ~ 1,
                    adj.P.Val_control < 0.05 & logFC_control < 0 ~ -1,
                    .default = 0
                ),
                outcome_heparin = case_when(
                    adj.P.Val_heparin < 0.05 & logFC_heparin > 0 ~ 1,
                    adj.P.Val_heparin < 0.05 & logFC_heparin < 0 ~ -1,
                    .default = 0
                )
            ) %>%
            mutate(
                outcome_combined = case_when(
                    outcome_control != 0 & outcome_heparin == 0 ~ "control",
                    outcome_control == 0 & outcome_heparin != 0 ~ "heparin",
                    outcome_control != 0 & outcome_heparin != 0 ~ "both",
                    outcome_control == 0 & outcome_heparin == 0 ~ "none"
                )
            ) %>%
            mutate(
                outcome_combined = outcome_combined %>%
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
            )
    }
)

list_plots <- map2(
    list_results,
    names_conditions_interest,
    \(results, names) {
        # Subset results to gene sets
        plots <- results %>%
            mutate(symbol_repel = case_when(
                outcome_combined %in% c("Both", "ns") ~ NA,
                .default = set
            )) %>%
            ggplot(
                aes(
                    x = logFC_control,
                    y = logFC_heparin,
                    fill = outcome_combined,
                    label = symbol_repel
                )
            ) +
            geom_point(
                shape = "circle filled",
                alpha = 0.6,
                stroke = 0.5,
                size = 2
            ) +
            geom_text_repel(
                max.overlaps = 50,
                force_pull = 0.5,
                min.segment.length = 0
            ) +
            geom_quadrant_lines(linetype = "dotted") +
            scale_fill_manual(values = c(
                "Control only" = "blue",
                "Heparin only" = "red",
                "Both" = "purple",
                "ns" = "grey"
            )) +
            scale_x_continuous(limits = symmetric_limits(c(-2.6, 2.6))) +
            scale_y_continuous(limits = symmetric_limits(c(-2.6, 2.6))) +
            coord_obs_pred() +
            theme_classic() +
            labs(
                x = str_c("deltaES ", names[[1]]),
                y = str_c("deltaES ", names[[2]]),
                fill = "Outcome"
            )

        # Return data
        return(plots)
    }
)
