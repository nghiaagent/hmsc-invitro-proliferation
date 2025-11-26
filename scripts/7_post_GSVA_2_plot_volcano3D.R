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

# Define outcome variable
outcome <- phenoData(quant_gsva[[1]])$condition_ID %>%
    factor(
        levels = c(
            "P5D3Untreated",
            "P7D3Untreated",
            "P13D3Untreated"
        ),
        labels = c(
            "5",
            "7",
            "13"
        )
    )

# Define axis specs
breaks <- seq(
    from = 0,
    to = 3.5,
    by = 0.5
)

# Derive GSVA results

## Expression matrix
data_gsva <- map(
    quant_gsva,
    \(eset) {
        exprs(eset)
    }
)

## Supply pvalues and padj
polar_pvals_gsva <- map(
    fit_gsva,
    \(fit_contrasts) {
        cbind(
            topTable(fit_contrasts, number = Inf, sort.by = "none")$P.Value,
            topTable(fit_contrasts, number = Inf, sort.by = "none", coef = 13)$P.Value,
            topTable(fit_contrasts, number = Inf, sort.by = "none", coef = 15)$P.Value,
            topTable(fit_contrasts, number = Inf, sort.by = "none", coef = 14)$P.Value
        )
    }
)

polar_padj_gsva <- map(
    fit_gsva,
    \(fit_contrasts) {
        cbind(
            topTable(fit_contrasts, number = Inf, sort.by = "none")$adj.P.Val,
            topTable(fit_contrasts, number = Inf, sort.by = "none", coef = 13)$adj.P.Val,
            topTable(fit_contrasts, number = Inf, sort.by = "none", coef = 15)$adj.P.Val,
            topTable(fit_contrasts, number = Inf, sort.by = "none", coef = 14)$adj.P.Val
        )
    }
)

# Construct volcano3d object
polar_manual_gsva <- pmap(
    list(
        data_gsva,
        polar_pvals_gsva,
        polar_padj_gsva
    ),
    \(data_gsva, polar_pvals_gsva, polar_padj_gsva) {
        polar_coords(
            outcome = outcome,
            data = t(data_gsva),
            pvals = polar_pvals_gsva,
            padj = polar_padj_gsva,
            scheme = palette_volcano3d
        )
    }
)

# Add rownames to padj and pvals
polar_manual_gsva <- map2(
    polar_manual_gsva,
    fit_gsva,
    \(polar_manual, fit_contrasts) {
        rownames(polar_manual@pvals) <- rownames(
            topTable(
                fit_contrasts,
                number = Inf,
                sort.by = "none"
            )
        )
        rownames(polar_manual@padj) <- rownames(
            topTable(
                fit_contrasts,
                number = Inf,
                sort.by = "none"
            )
        )

        return(polar_manual)
    }
)

# Plot
volcano3d <- map(
    polar_manual_gsva,
    \(polar_manual) {
        volcano3D(
            polar_manual,
            type = 1,
            label_size = 24,
            z_axis_title_size = 20,
            radial_axis_title_size = 20
        )
    }
)

radial_plotly <- map(
    polar_manual_gsva,
    \(polar_manual) {
        radial_plotly(
            polar_manual,
            type = 1,
            r_axis_ticks = breaks
        )
    }
)

radial_ggplot <- map(
    polar_manual_gsva,
    \(polar_manual) {
        radial_ggplot(
            polar_manual,
            type = 1,
            r_axis_ticks = breaks
        )
    }
)

# Save data
saveRDS(
    polar_manual_gsva,
    file = here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "polar_manual_gsva.RDS"
    )
)

saveRDS(
    list(
        volcano3d,
        radial_plotly,
        radial_ggplot
    ),
    file = here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "volcano3d_gsva_passages_day3.RDS"
    )
)
