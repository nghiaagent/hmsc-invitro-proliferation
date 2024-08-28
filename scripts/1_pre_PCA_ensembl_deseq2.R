quant_big <- readRDS("./output/data_expression/pre_DGE/quant_cDNA_deseq.RDS")

# Try RUVg - normalise by endogenous controls

## Define endo controls

list_genes <- rownames(quant_big)

list_endogenous <- c("ACTB",
                     "B2M",
                     "EEF1A1",
                     "GAPDH",
                     "RPLP0",
                     "RPL13A",
                     "TBP",
                     "YWHAZ") %>%
  AnnotationDbi::select(org.Hs.eg.db,
                        keys = .,
                        keytype = "SYMBOL",
                        columns = "ENSEMBL") %>%
  .$ENSEMBL

## Convert object to SeqExpressionSet, run RUVg

quant_RUVg <- newSeqExpressionSet(counts = counts(quant_big),
                                  phenoData = colData(quant_big) %>% as.data.frame()) %>%
  RUVg(., list_endogenous, k = 4)

## Convert back to DESeqDataSet for plotting

quant_vsd_RUVg <- quant_big
x <- normCounts(quant_RUVg)
mode(x) <- "integer"
counts(quant_vsd_RUVg) <- x
rm(x)
quant_vsd_RUVg %<>% vst(., blind = FALSE)

# Try RUVs - normalise by replicates

## Define replicates

groups <- str_c(as.data.frame(colData(quant_big))$cell_line,
                as.data.frame(colData(quant_big))$condition_ID) %>%
  makeGroups()

## Convert object to SeqExpressionSet, run RUVs

quant_RUVs <- newSeqExpressionSet(counts = counts(quant_big),
                                  phenoData = colData(quant_big) %>% as.data.frame()) %>%
  RUVs(., list_genes, k = 4, groups)

## Convert back to DESeqDataSet for plotting

quant_vsd_RUVs <- quant_big
x <- normCounts(quant_RUVs)
mode(x) <- "integer"
counts(quant_vsd_RUVs) <- x
rm(x)
quant_vsd_RUVs %<>% vst(., blind = FALSE)

# Plot PCA

plotPCA_paired <- function(object,
                           intgroup = "condition",
                           ntop = 500,
                           returnData = FALSE,
                           pcsToUse = 1:2) {
  message(paste0("using ntop=", ntop, " top features by variance"))
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev ^ 2 / sum(pca$sdev ^ 2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  pcs <- paste0("PC", pcsToUse)
  d <- data.frame(
    V1 = pca$x[, pcsToUse[1]],
    V2 = pca$x[, pcsToUse[2]],
    group = group,
    intgroup.df,
    name = colnames(object)
  )
  colnames(d)[1:2] <- pcs
  if (returnData) {
    attr(d, "percentVar") <- percentVar[pcsToUse]
    return(d)
  }
  ggplot(data = d, aes_string(x = pcs[1], y = pcs[2], color = "group")) +
    geom_point(size = 3) +
    xlab(paste0(pcs[1], ": ", round(percentVar[pcsToUse[1]] * 100), "% variance")) +
    ylab(paste0(pcs[2], ": ", round(percentVar[pcsToUse[2]] * 100), "% variance")) +
    coord_fixed() +
    scale_color_brewer(palette = "Paired")
}

plotPCA_all <- function(intgroup) {
  plot_grid(
    plotPCA_paired(vst(quant_big, blind = FALSE), intgroup = intgroup),
    plotPCA_paired(quant_vsd_RUVg, intgroup = intgroup),
    plotPCA_paired(quant_vsd_RUVs, intgroup = intgroup)
  )
}

quant_big_batchcor <- quant_big

counts(quant_big_batchcor) <- sva::ComBat_seq(counts(quant_big_batchcor),
                                              batch = colData(quant_big_batchcor)$run_date
                                              # ,
                                              # covar_mod = model.matrix(~ Passage, data = table_samples)
                                              ) %>%
  `storage.mode<-`(., "integer")

cor <- vst(quant_big_batchcor, blind = FALSE) %$%
  plot_grid(
    plotPCA_paired(., intgroup = "run_date"),
    plotPCA_paired(., intgroup = "cell_line"),
    plotPCA_paired(., intgroup = "condition_ID"),
    plotPCA_paired(., intgroup = "Passage"),
    plotPCA_paired(., intgroup = "Day"),
    plotPCA_paired(., intgroup = "Treatment"),
    align = "v"
  )

uncor <- vst(quant_big, blind = FALSE) %$%
  plot_grid(
    plotPCA_paired(., intgroup = "run_date"),
    plotPCA_paired(., intgroup = "cell_line"),
    plotPCA_paired(., intgroup = "condition_ID"),
    plotPCA_paired(., intgroup = "Passage"),
    plotPCA_paired(., intgroup = "Day"),
    plotPCA_paired(., intgroup = "Treatment"),
    align = "v"
  )