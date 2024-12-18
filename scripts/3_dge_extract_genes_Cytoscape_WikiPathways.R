# Prepare dataset. Run DGE first.
# Swap ENSEMBL ID with ENTREZ ID
# If multiple ENS genes share one ENTREZ ID, remove gene with lower average expression
# Remove genes with no Entrez ID

fit_contrasts <-
    readRDS(file = "./output/data_expression/post_DGE/fit_contrasts.RDS")

order <- order(fit_contrasts$Amean, decreasing = TRUE)

# Extract logFC between untreated cells

table_logFC <- topTable(
    fit_contrasts,
    coef = c(13, 14, 15),
    sort.by = "none",
    n = Inf
)

table_logFC <- select(
    table_logFC,
    c(
        GENEID,
        GENENAME,
        P7vsP5_UT_D3,
        P13vsP7_UT_D3,
        P13vsP5_UT_D3
    )
)

write_csv(table_logFC,
    file = "./output/list_genes/logFC_betweenpassages_Cytoscape.csv"
)
