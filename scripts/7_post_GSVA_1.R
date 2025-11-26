# Load data

rlog_deseq2 <- readRDS(
    file = here::here(
        "output",
        "data_expression",
        "post_DGE",
        "rlog_deseq2.RDS"
    )
)

# Convert to ESet, convert identifiers to ENTREZID
# Sort dataset by decreasing expression
# Remove genes with duplicate ENTREZ ID

order <- order(
    rowRanges(rlog_deseq2)$baseMean,
    decreasing = TRUE
)

rlog_deseq2 <- rlog_deseq2[order, ]

rlog_deseq2 <- rlog_deseq2[
    rowRanges(rlog_deseq2) %>%
        as.data.frame() %$%
        map(
            .$entrezid,
            \(x) !is.na(x)[[1]]
        ) %>%
        unlist(),
]

idx <- rowRanges(rlog_deseq2) %>%
    as.data.frame() %$%
    map(
        .$entrezid,
        \(x) x[[1]]
    ) %>%
    unlist()

rlog_deseq2 <- rlog_deseq2[!duplicated(idx), ]
idx <- idx[!duplicated(idx)]

# Rename rownames to entrezid

rlog_deseq2_counts <- assay(rlog_deseq2)
rlog_deseq2_rowranges <- rowRanges(rlog_deseq2) %>%
    as.data.frame()
rlog_deseq2_coldata <- colData(rlog_deseq2) %>%
    as.data.frame()

rownames(rlog_deseq2_counts) <- idx
rownames(rlog_deseq2_rowranges) <- idx

# Build eset object

quant_eset <- ExpressionSet(
    assayData = rlog_deseq2_counts,
    phenoData = AnnotatedDataFrame(rlog_deseq2_coldata),
    featureData = AnnotatedDataFrame(rlog_deseq2_rowranges)
)

# Define design and contrasts for GSVA + limma
design <- model.matrix(~ condition_ID + cell_line,
    data = rlog_deseq2_coldata
)

colnames(design) <- make.names(colnames(design))

matrix_contrasts <- makeContrasts(

    ## Coefs 1 - 6: Treatment at each timepoint
    Trt_P5_D3 = condition_IDP5D3Treated - 0,
    Trt_P5_D5 = condition_IDP5D5Treated - condition_IDP5D5Untreated,
    Trt_P7_D3 = condition_IDP7D3Treated - condition_IDP7D3Untreated,
    Trt_P7_D5 = condition_IDP7D5Treated - condition_IDP7D5Untreated,
    Trt_P13_D3 = condition_IDP13D3Treated - condition_IDP13D3Untreated,
    Trt_P13_D5 = condition_IDP13D5Treated - condition_IDP13D5Untreated,

    ## Coefs 7 - 12: Day at each timepoint x treatment
    D5vsD3_UT_P5 = condition_IDP5D5Untreated - 0,
    D5vsD3_UT_P7 = condition_IDP7D5Untreated - condition_IDP7D3Untreated,
    D5vsD3_UT_P13 = condition_IDP13D5Untreated - condition_IDP13D3Untreated,
    D5vsD3_T_P5 = condition_IDP5D5Treated - condition_IDP5D3Treated,
    D5vsD3_T_P7 = condition_IDP7D5Treated - condition_IDP7D3Treated,
    D5vsD3_T_P13 = condition_IDP13D5Treated - condition_IDP13D3Treated,

    ## Coefs 13 - 24: Passage at each day x treatment
    P7vsP5_UT_D3 = condition_IDP7D3Untreated - 0,
    P13vsP7_UT_D3 = condition_IDP13D3Untreated - condition_IDP7D3Untreated,
    P13vsP5_UT_D3 = condition_IDP13D3Untreated - 0,
    P7vsP5_T_D3 = condition_IDP7D3Treated - condition_IDP5D3Treated,
    P13vsP7_T_D3 = condition_IDP13D3Treated - condition_IDP7D3Treated,
    P13vsP5_T_D3 = condition_IDP13D3Treated - condition_IDP5D3Treated,
    P7vsP5_UT_D5 = condition_IDP7D5Untreated - condition_IDP5D5Untreated,
    P13vsP7_UT_D5 = condition_IDP13D5Untreated - condition_IDP7D5Untreated,
    P13vsP5_UT_D5 = condition_IDP13D5Untreated - condition_IDP5D5Untreated,
    P7vsP5_T_D5 = condition_IDP7D5Treated - condition_IDP5D5Treated,
    P13vsP7_T_D5 = condition_IDP13D5Treated - condition_IDP7D5Treated,
    P13vsP5_T_D5 = condition_IDP13D5Treated - condition_IDP5D5Treated,
    levels = design
)

# Perform GSVA
## Make GSVA NES object

quant_gsva_db <- map(
    list_gmt[1:6],
    \(x) {
        gsvaParam(
            quant_eset,
            x,
            minSize = 5,
            maxSize = 500,
            kcdf = "Gaussian"
        ) %>%
            gsva()
    },
    .progress = TRUE
)

quant_gsva_wgcna <- map(
    list_gmt[7],
    \(x) {
        gsvaParam(
            quant_eset,
            x,
            kcdf = "Gaussian"
        ) %>%
            gsva()
    },
    .progress = TRUE
)

quant_gsva <- c(
    quant_gsva_db,
    quant_gsva_wgcna
)

## Perform model fit
fit_gsva <- map(
    quant_gsva,
    \(x) {
        lmFit(x, design) %>%
            eBayes() %>%
            contrasts.fit(matrix_contrasts) %>%
            eBayes()
    },
    .progress = TRUE
)

# Save data

saveRDS(
    fit_gsva,
    here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "GSVA_results.RDS"
    )
)

saveRDS(
    quant_gsva,
    here::here(
        "output",
        "data_enrichment",
        "GSVA",
        "quant_GSVA.RDS"
    )
)
