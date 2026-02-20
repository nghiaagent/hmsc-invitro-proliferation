# Load data

## GCN with all samples
gcn <- readRDS(
    file = here::here(
        "output",
        "data_WGCNA",
        "WGCNA_allsamples",
        "GCN_output.RDS"
    )
)

## Add module assignment to gene list
gcn$genes$module <- gcn$net$colors %>%
    labels2colors()

gcn$genes <- gcn$genes %>%
    rownames_to_column(var = "entrezid_unique")

## Extract gene lists
## entrez for GSVA
## Symbol for String/PPI

gcn_genelists_entrez <- split(
    gcn$genes,
    gcn$genes$module,
    drop = TRUE
) %>%
    lapply(
        function(x) {
            select(x, entrezid_unique) %>%
                as_vector() %>%
                unname()
        }
    )

gcn_genelists_symbol <- split(
    gcn$genes,
    gcn$genes$module,
    drop = TRUE
) %>%
    lapply(
        function(x) {
            select(x, symbol) %>%
                as_vector() %>%
                unname()
        }
    )

names(gcn_genelists_entrez) <- str_c(
    "WGCNA_allsamples_",
    names(gcn_genelists_entrez)
)

names(gcn_genelists_symbol) <- str_c(
    "WGCNA_allsamples_",
    names(gcn_genelists_symbol)
)

# Define function to write list to GMT file
### Function by Levi Waldron.
### R list object that will be converted to GMT file.  Each element
### should contain a vector of gene names, and the names of the
### elements will used for the gene set names
### Output file name for .gmt file

write_gmt <- function(object, fname) {
    if (class(object) != "list") {
        stop("object should be of class 'list'")
    }

    if (file.exists(fname)) {
        unlink(fname)
    }

    for (iElement in 1:length(object)) { # nolint
        write.table(
            t(
                c(
                    make.names(
                        rep(
                            names(object)[iElement],
                            2
                        )
                    ),
                    object[[iElement]]
                )
            ),
            sep = "\t",
            quote = FALSE,
            file = fname,
            append = TRUE,
            col.names = FALSE,
            row.names = FALSE
        )
    }
}

# Write gene list in ENTREZID to GMT file

write_gmt(
    gcn_genelists_entrez,
    here::here(
        ".",
        "input",
        "genesets",
        "gcn_sets_WGCNA.gmt"
    )
)

# Write gene list with GENENAME to text files for STRING

imap(
    list(
        "turquoise_genes" = "WGCNA_allsamples_turquoise",
        "pink_genes" = "WGCNA_allsamples_pink",
        "blue_genes" = "WGCNA_allsamples_blue",
        "red_genes" = "WGCNA_allsamples_red",
        "cyan_genes" = "WGCNA_allsamples_cyan"
    ),
    \(table, name) {
        fwrite(
            list(gcn_genelists_symbol[[table]]),
            file = here::here(
                "output",
                "data_WGCNA",
                "WGCNA_allsamples",
                str_c(name, ".txt")
            ),
            na = ""
        )
    }
)
