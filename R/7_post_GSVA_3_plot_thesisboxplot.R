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

#
