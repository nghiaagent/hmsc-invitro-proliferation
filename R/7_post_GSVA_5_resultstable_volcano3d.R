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

polar_manual_gsva <- readRDS(
    file = here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "polar_manual_gsva.RDS"
    )
)

# Define table specs
ntop <- 5

# Define coefficients for extracting genes
coef_volcano3d <- list(
    up = list(
        p5  = "5+",
        p7  = "7+",
        p13 = "1+"
    ),
    down = list(
        p5  = "7+1+",
        p7  = "5+1+",
        p13 = "5+7+"
    )
)

coef_toptable <- list(
    p5  = c(p5_p7 = "P7vsP5_UT_D3", p5_p13 = "P13vsP5_UT_D3"),
    p7  = c(p5_p7 = "P7vsP5_UT_D3", p7_p13 = "P13vsP7_UT_D3"),
    p13 = c(p7_p13 = "P13vsP7_UT_D3", p5_p13 = "P13vsP5_UT_D3")
)

# Define collections to consider

collections_interest <- c(
    "h",
    "c2_cgp",
    "c2_cp",
    "GOBP",
    "GOCC",
    "GOMF"
)

# Extract gene set IDs relevant
genesetid_gsva <- map(
    polar_manual_gsva[1:6],
    \(polar_manual) {
        map(
            coef_volcano3d,
            \(x) {
                map(
                    x,
                    \(coefs) {
                        significance_subset(
                            polar_manual,
                            significance = coefs,
                            output = "pvals"
                        ) %>%
                            rownames()
                    }
                )
            }
        )
    }
)

genesetid_wgcna <- map(
    coef_volcano3d,
    \(x) {
        map(
            x,
            \(coefs) {
                polar_manual_gsva[["WGCNA"]]@df$scaled %>%
                    filter(lab == coefs) %>%
                    rownames()
            }
        )
    }
)

genesetid_gsva <- c(genesetid_gsva, list("WGCNA" = genesetid_wgcna))

# Merge fit_gsva lists to get list of
# genesets up/downregulated at each passage compared to the rest
toptable_paired <- imap(
    fit_gsva,
    \(fit, name) {
        map(
            coef_toptable,
            \(contrast) {
                extract_joined_results_limma(
                    fit,
                    contrast_1 = contrast[[1]],
                    contrast_2 = contrast[[2]]
                ) %>%
                    mutate(collection = name)
            }
        )
    }
)

# Select for gene sets highlighted in 3D volcano plots
toptable_ordered <- map2(
    toptable_paired,
    genesetid_gsva,
    \(toptable_collection, genesetid_collection) {
        map(
            genesetid_collection,
            \(genesetid) {
                map2(
                    toptable_collection,
                    genesetid,
                    \(top, genesetid) {
                        filter(top, geneset %in% genesetid) %>%
                            dplyr::arrange(., desc(.[[2]]))
                    }
                )
            }
        )
    }
)

# Merge into one table for export to thesis
toptable_out <- toptable_ordered %$%
    list(
        p5 = rbind(
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["down"]][["p5"]],
                        P7vsP5_UT_D3 > 0,
                        P13vsP5_UT_D3 > 0
                    ) %>%
                        drop_na() %>%
                        slice_head(n = ntop)
                }
            )),
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["up"]][["p5"]],
                        P7vsP5_UT_D3 < 0,
                        P13vsP5_UT_D3 < 0
                    ) %>%
                        drop_na() %>%
                        slice_tail(n = ntop)
                }
            ))
        ),
        p7 = rbind(
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["up"]][["p7"]],
                        P7vsP5_UT_D3 > 0,
                        P13vsP7_UT_D3 < 0
                    ) %>%
                        drop_na() %>%
                        slice_head(n = ntop)
                }
            )),
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["down"]][["p7"]],
                        P7vsP5_UT_D3 < 0,
                        P13vsP7_UT_D3 > 0
                    ) %>%
                        drop_na() %>%
                        slice_tail(n = ntop)
                }
            ))
        ),
        p13 = rbind(
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["up"]][["p13"]],
                        P13vsP7_UT_D3 > 0,
                        P13vsP5_UT_D3 > 0
                    ) %>%
                        drop_na() %>%
                        slice_head(n = ntop)
                }
            )),
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["down"]][["p13"]],
                        P13vsP7_UT_D3 < 0,
                        P13vsP5_UT_D3 < 0
                    ) %>%
                        drop_na() %>%
                        slice_tail(n = ntop)
                }
            ))
        )
    )

toptable_out_big <- toptable_ordered %$%
    list(
        p5 = rbind(
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["down"]][["p5"]],
                        P7vsP5_UT_D3 > 0,
                        P13vsP5_UT_D3 > 0
                    ) %>%
                        drop_na()
                }
            )),
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["up"]][["p5"]],
                        P7vsP5_UT_D3 < 0,
                        P13vsP5_UT_D3 < 0
                    ) %>%
                        drop_na()
                }
            ))
        ),
        p7 = rbind(
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["up"]][["p7"]],
                        P7vsP5_UT_D3 > 0,
                        P13vsP7_UT_D3 < 0
                    ) %>%
                        drop_na()
                }
            )),
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["down"]][["p7"]],
                        P7vsP5_UT_D3 < 0,
                        P13vsP7_UT_D3 > 0
                    ) %>%
                        drop_na()
                }
            ))
        ),
        p13 = rbind(
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["up"]][["p13"]],
                        P13vsP7_UT_D3 > 0,
                        P13vsP5_UT_D3 > 0
                    ) %>%
                        drop_na()
                }
            )),
            bind_rows(map(
                collections_interest,
                \(collection) {
                    dplyr::filter(
                        .[[collection]][["down"]][["p13"]],
                        P13vsP7_UT_D3 < 0,
                        P13vsP5_UT_D3 < 0
                    ) %>%
                        drop_na()
                }
            ))
        )
    )


# Export data
write.xlsx(
    x = toptable_out,
    asTable = TRUE,
    file <- here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "toptable_volcano3d.xlsx"
    )
)

write.xlsx(
    x = toptable_out_big,
    asTable = TRUE,
    file <- here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "toptable_volcano3d_big.xlsx"
    )
)
