# Load data
# Transpose logCPM matrix for WGCNA
# Convert identifiers from ENSEMBL to ENTREZ (to help with gene set testing)
# Remove genes with duplicate ENTREZID

quant_deseq2 <- readRDS(
    file = here::here(
        "output",
        "data_expression",
        "post_DGE",
        "quant_deseq2_batchcor.RDS"
    )
)

rlog_deseq2 <- readRDS(
    file = here::here(
        "output",
        "data_expression",
        "post_DGE",
        "rlog_deseq2.RDS"
    )
)

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

# Create matrix for GCN

gcn_rlog <- t(rlog_deseq2_counts) %>%
    as.data.frame()

# Check for low quality genes and samples
## Should be all OK
## If not, remove offending genes and samples from data

gcn_qual <- goodSamplesGenes(gcn_rlog)

if (!gcn_qual$allOK) {
    # Print genes that are removed
    if (sum(!gcn_qual$goodGenes) > 0) {
        printFlush(
            paste(
                "Removing genes:",
                paste(
                    names(datExpr0)[!gcn_qual$goodGenes],
                    collapse = ", "
                )
            )
        )
    }

    # Print samples that are removed
    if (sum(!gcn_qual$goodSamples) > 0) {
        printFlush(
            paste(
                "Removing samples:",
                paste(
                    rownames(datExpr0)[!gcn_qual$goodSamples],
                    collapse = ", "
                )
            )
        )
    }

    # Remove the offending genes and samples from the data:
    gcn_rlog <- gcn_rlog[gcn_qual$goodSamples, gcn_qual$goodGenes]
}

# Perform sample clustering to detect outliers
## One outlier detected and removed (hMSC-21558_P13_D5_10-1)

metadata_select <- c(
    "cell_line",
    "Passage",
    "Day",
    "Treatment",
    "run_date"
)

png(
    here::here(
        "output",
        "plots_WGCNA",
        "WGCNA_allsamples",
        "Sample clustering.png"
    ),
    res = 150,
    width = 40,
    height = 25,
    units = "cm"
)

plotDendroAndColors(
    dendro = hclust(
        dist(gcn_rlog),
        method = "average"
    ),
    colors = labels2colors(dplyr::select(
        rlog_deseq2_coldata,
        all_of(metadata_select)
    )),
    groupLabels = names(dplyr::select(
        rlog_deseq2_coldata,
        all_of(metadata_select)
    )),
    main = "Hierarchical clustering of samples"
)

dev.off()

## Remove detected outlier from gcn_rlog and rlog_deseq2_coldata

samples_exclude <- c(
    "hMSC-21558_P13_D5_10-1"
)

gcn_rlog <- gcn_rlog[!rownames(gcn_rlog) %in% samples_exclude, ]
rlog_deseq2_coldata <- rlog_deseq2_coldata %>%
    dplyr::filter(
        !names %in% samples_exclude
    )

## Hierarchical clustering after removal of outliers

png(
    here::here(
        "output",
        "plots_WGCNA",
        "WGCNA_allsamples",
        "Sample clustering (after outlier removal).png"
    ),
    res = 150,
    width = 40,
    height = 25,
    units = "cm"
)

plotDendroAndColors(
    dendro = hclust(
        dist(gcn_rlog),
        method = "average"
    ),
    colors = labels2colors(dplyr::select(
        rlog_deseq2_coldata,
        all_of(metadata_select)
    )),
    groupLabels = names(dplyr::select(
        rlog_deseq2_coldata,
        all_of(metadata_select)
    )),
    main = "Hierarchical clustering of samples"
)

dev.off()

# Select soft thresholding power
## Calculate scale-free topology fit index

powers <- c(
    c(1:10),
    seq(
        from = 12,
        to = 20,
        by = 2
    )
)

gcn_sft <- pickSoftThreshold(
    gcn_rlog,
    powerVector = powers,
    verbose = 5
)

## Plot scale-free topology fit index and mean connectivity

png(
    here::here(
        "output",
        "plots_WGCNA",
        "WGCNA_allsamples",
        "Soft thresholding index selection.png"
    ),
    res = 300,
    width = 20,
    height = 10,
    units = "cm"
)

par(mfrow = c(1, 2))
cex1 <- 0.9

### SFT index

plot(
    gcn_sft$fitIndices[, 1],
    -sign(gcn_sft$fitIndices[, 3]) * gcn_sft$fitIndices[, 2],
    xlab = "Soft Threshold (power)",
    ylab = "Scale Free Topology Model Fit,signed R^2",
    type = "n",
    main = paste("Scale independence")
)

text(
    gcn_sft$fitIndices[, 1],
    -sign(gcn_sft$fitIndices[, 3]) * gcn_sft$fitIndices[, 2],
    labels = powers,
    cex = cex1,
    col = "red"
)

abline(h = 0.90, col = "red")

### Mean connectivity

plot(
    gcn_sft$fitIndices[, 1],
    gcn_sft$fitIndices[, 5],
    xlab = "Soft Threshold (power)",
    ylab = "Mean Connectivity",
    type = "n",
    main = paste("Mean connectivity")
)

text(
    gcn_sft$fitIndices[, 1],
    gcn_sft$fitIndices[, 5],
    labels = powers,
    cex = cex1,
    col = "red"
)

dev.off()

# Construct GCN object

gcn <- list(
    E = gcn_rlog,
    targets = rlog_deseq2_coldata %>%
        mutate(
            cell_line = as.numeric(cell_line),
            Passage = as.numeric(Passage),
            Day = as.numeric(Day),
            Treatment = as.numeric(Treatment),
            run_date = as.numeric(factor(run_date)),
            timepoint_ID = as.numeric(timepoint_ID)
        ) %>%
        select(
            cell_line,
            Passage,
            Day,
            Treatment,
            run_date,
            timepoint_ID
        ),
    genes = rlog_deseq2_rowranges,
    sft = gcn_sft
)

saveRDS(gcn,
    file = here::here(
        "output",
        "data_WGCNA",
        "WGCNA_allsamples",
        "GCN_input.Rdata"
    )
)
