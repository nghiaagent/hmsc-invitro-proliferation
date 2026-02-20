# Extract gene names from DESeq fit for ORA using g:Profiler
# All contrasts

# Load data
quant_deseq2 <- readRDS(
    file = here::here(
        "output", "data_expression", "post_DGE",
        "quant_deseq2_batchcor.RDS"
    )
)

results_lfcshrink <- readRDS(
    file = here::here(
        "output",
        "data_expression",
        "post_DGE",
        "results_deseq2_lfcshrink.RDS"
    )
)

# Format table object, get only gene ENSEMBL IDs, ENTREZ IDs, and symbol
# Enter ENSEMBL IDs or ENTREZ IDs into g:Profiler
# Enter symbol into STRING
results_lfcshrink_format <- map(
    results_lfcshrink,
    \(x) {
        extract_topgenes(
            results = x,
            dds = quant_deseq2,
            ntop = Inf,
            signif_only = TRUE
        )
    },
    .progress = TRUE
)

results_ensembl_ids <- map(
    results_lfcshrink_format,
    \(x) x$`ENSEMBL ID`
)

results_entrez_ids <- map(
    results_lfcshrink_format,
    \(x) x$`ENTREZ ID`
)

results_symbol <- map(
    results_lfcshrink_format,
    \(x) x$Symbol
)

# Export data

## Export ENSEMBL IDs
imap(
    results_ensembl_ids,
    \(x, idx) {
        fwrite(
            list(x),
            file = here::here(
                "output",
                "list_genes",
                str_c(idx, "_ENSEMBL.txt")
            )
        )
    }
)

## Export gene names
imap(
    results_symbol,
    \(x, idx) {
        fwrite(
            list(x),
            file = here::here(
                "output",
                "list_genes",
                str_c(idx, "_symbol.txt")
            )
        )
    }
)

## Export ENTREZ ID

imap(
    results_entrez_ids,
    \(x, idx) {
        fwrite(
            list(x),
            file = here::here(
                "output",
                "list_genes",
                str_c(idx, "_ENTREZ.txt")
            )
        )
    }
)
