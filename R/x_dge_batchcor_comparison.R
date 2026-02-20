# Compare F-statistics produced by
# Non batch corrected data, with modelled batch (fit_contrasts_1)
# Batch corrected data, without modelling for batch (fit_contrasts_2)
# Batch corrected data, with modelling for batch (fit_contrasts_3)
# Run DGE script first

design_nobatch <- model.matrix(~ condition_ID + cell_line,
    data = table_design
)

fit_contrasts_1 <-
    voom(quant_DGE_clean,
        design,
        plot = TRUE
    ) %>%
    lmFit(design) %>%
    eBayes() %>%
    contrasts.fit(matrix_contrasts) %>%
    eBayes()

fit_contrasts_2 <-
    voom(quant_DGE_clean_batchcor,
        design_nobatch,
        plot = TRUE
    ) %>%
    lmFit(design) %>%
    eBayes() %>%
    contrasts.fit(matrix_contrasts) %>%
    eBayes()

fit_contrasts_3 <- fit_contrasts

pvals_1 <- topTable(fit_contrasts_1, coef = 15, number = Inf, sort.by = "none") %>%
    select(GENEID, P.Value)

pvals_2 <- topTable(fit_contrasts_2, coef = 15, number = Inf, sort.by = "none") %>%
    select(GENEID, P.Value)

pvals_3 <- topTable(fit_contrasts_3, coef = 15, number = Inf, sort.by = "none") %>%
    select(GENEID, P.Value)

pvals_plot <- left_join(pvals_1,
    pvals_3,
    by = join_by(GENEID == GENEID),
    suffix = c("_1", "_2")
) %>%
    pivot_longer(
        cols = !GENEID,
        names_to = "type",
        names_prefix = "P.Value_"
    ) %>%
    mutate(type = factor(type,
        levels = c("1", "2"),
        labels = c(
            "Model only",
            "Model and correct"
        )
    ))

plot <-
    ggplot(
        data = pvals_plot,
        aes(
            x = value,
            colors = type,
            fill = type
        )
    ) +
    geom_density(alpha = 0.2)

plot_grid(
    plot +
        facet_wrap(~type) +
        theme(legend.position = "none"),
    plot,
    nrow = 1
) %>%
    ggsave(
        filename = file.path(
            "output",
            "plots_QC",
            "F-statistics between batch correction methods.png"
        ),
        plot = .,
        scale = 0.8,
        width = 12,
        height = 6
    )
