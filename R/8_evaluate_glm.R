here::i_am("R/8_evaluate_glm.R")

########################
# Evaluate GLM methods (edgeR vs limma vs DESeq2)
########################

# Import packages
library(colorspace)
library(conflicted)
library(DESeq2)
library(edgeR)
library(ggVennDiagram)
library(here)
library(limma)
library(magrittr)
library(SummarizedExperiment)
library(tidyverse)

# Set variables
## Cutoffs for volcano and quadrant plots
cutoff_logfc <- 6
cutoff_logfc_neg <- cutoff_logfc * -1
cutoff_padj <- 1e-15

## Palette for quadrant plots
palette_quadrant_evaluate_glm <- list(
  deseq2_limma = c(
    "limma only" = palette()[[4]],
    "DESeq2 only" = palette()[[2]],
    "Both" = palette()[[6]],
    "None" = palette()[[8]] %>% darken(0.4)
  ),
  edger_limma = c(
    "limma only" = palette()[[4]],
    "edgeR only" = palette()[[2]],
    "Both" = palette()[[6]],
    "None" = palette()[[8]] %>% darken(0.4)
  ),
  deseq2_edger = c(
    "edgeR only" = palette()[[4]],
    "DESeq2 only" = palette()[[2]],
    "Both" = palette()[[6]],
    "None" = palette()[[8]] %>% darken(0.4)
  )
)

# Load data
## DESeq2
quant_deseq2 <- readRDS(here::here(
  "output",
  "data_expression",
  "pre_DGE",
  "quant_cDNA_deseq.RDS"
))

## Select contrasts for evaluation
list_contrasts_deseq2_pilot <- list_contrasts_deseq2[
  c(
    "P7vsP5_UT_D3",
    "P13vsP7_UT_D3",
    "P13vsP5_UT_D3"
  )
]

## Get gene names
gene_names <- SummarizedExperiment::rowRanges(quant_deseq2) %>%
  as.data.frame() %>%
  dplyr::select(
    gene_id,
    gene_name
  )

## edgeR
### Convert quant_deseq2 to DGEList
quant_edger <- quant_deseq2 %$%
  edgeR::DGEList(
    counts = DESeq2::counts(.),
    samples = SummarizedExperiment::colData(.),
    group = SummarizedExperiment::colData(.)$condition_ID
  ) %>%
  edgeR::normLibSizes()

### Design matrix (also used for limma)
design <- model.matrix(
  ~ condition_ID + run_date + cell_line,
  data = quant_edger$samples
)

colnames(design) <- make.names(colnames(design))

### Contrast matrix (also used for limma)
contrasts_pilot <- limma::makeContrasts(
  ## Coefs 1 - 3: Passage at each day x treatment
  P7vsP5_UT_D3 = condition_IDP7D3Untreated - 0,
  P13vsP7_UT_D3 = condition_IDP13D3Untreated - condition_IDP7D3Untreated,
  P13vsP5_UT_D3 = condition_IDP13D3Untreated - 0,
  levels = design
)

list_coefs <- ncol(contrasts_pilot) %>%
  seq_len() %>%
  magrittr::set_names(colnames(contrasts_pilot))

## limma-voom
quant_limma <- quant_edger %>%
  limma::voom(design, plot = TRUE)

# Analyse with DESeq2
quant_deseq2 <- quant_deseq2 %>%
  DESeq2::DESeq()

## Obtain results
### No LFC shrinking
results_deseq2 <- list_contrasts_deseq2_pilot %>%
  purrr::map(\(contrast) {
    DESeq2::results(
      quant_deseq2,
      contrast = contrast,
      filterFun = ihw,
      alpha = 0.05
    )
  })

# Analyse with edgeR (glm)
quant_edger <- quant_edger %>%
  edgeR::estimateDisp(design, robust = TRUE)

fit_edger <- quant_edger %>%
  edgeR::glmQLFit(design, robust = TRUE)

results_edger <- list_coefs %>%
  purrr::map(\(coef) {
    edgeR::glmQLFTest(
      fit_edger,
      contrast = contrasts_pilot[, coef]
    )
  })

top_edger <- results_edger %>%
  purrr::map(\(x) {
    edgeR::topTags(
      x,
      n = Inf,
      sort.by = "none"
    )
  })

# Analyse with limma-voom
fit_limma_contrasts <- quant_limma %>%
  limma::lmFit(design) %>%
  limma::eBayes() %>%
  limma::contrasts.fit(contrasts_pilot) %>%
  limma::eBayes()

top_limma <- list_coefs %>%
  purrr::map(\(coef) {
    limma::topTable(
      fit_limma_contrasts,
      coef = coef,
      number = Inf,
      sort.by = "none"
    )
  })

# Get summary statistics
## DESeq2
results_deseq2 %>%
  purrr::map(\(x) summary(x))

## edgeR
list_coefs %>%
  purrr::map(\(coef) {
    edgeR::glmQLFTest(
      fit_edger,
      contrast = contrasts_pilot[, coef]
    ) %>%
      decideTests() %>%
      summary()
  })

## limma
fit_limma_contrasts %>%
  decideTests() %>%
  summary()

# Summarise results to one table per contrast
## Add gene column
top_edger <- top_edger %>%
  purrr::map(\(top) {
    top$table %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "ensembl_id") %>%
      tibble::as_tibble() %>%
      dplyr::rename_with(
        ~ stringr::str_c(.x, "_edger"),
        !ensembl_id
      )
  })

top_limma <- top_limma %>%
  purrr::map(\(top) {
    top %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "ensembl_id") %>%
      tibble::as_tibble() %>%
      dplyr::rename_with(
        ~ stringr::str_c(.x, "_limma"),
        !ensembl_id
      )
  })

results_deseq2 <- results_deseq2 %>%
  purrr::map(\(top) {
    top %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "ensembl_id") %>%
      tibble::as_tibble() %>%
      dplyr::rename_with(
        ~ stringr::str_c(.x, "_deseq2"),
        !ensembl_id
      )
  })

## Join into one table
top_joined <- results_deseq2 %>%
  purrr::map2(
    top_edger,
    \(deseq2, edger) {
      dplyr::full_join(
        deseq2,
        edger,
        by = dplyr::join_by("ensembl_id" == "ensembl_id")
      )
    }
  ) %>%
  # Join with limma data
  purrr::map2(
    top_limma,
    \(deseq2_edger, limma) {
      dplyr::full_join(
        deseq2_edger,
        limma,
        by = dplyr::join_by("ensembl_id" == "ensembl_id")
      )
    }
  ) %>%
  # Add gene name
  purrr::map(\(x) {
    dplyr::left_join(
      x,
      gene_names,
      by = dplyr::join_by("ensembl_id" == "gene_id")
    ) %>%
      dplyr::relocate(ensembl_id, gene_name)
  })

## Create table for volcano plot
top_joined_volcano <- top_joined %>%
  # Clip data
  purrr::map(\(results) {
    results <- results %>%
      # Clip logFC to threshold
      dplyr::mutate(
        log2FoldChange_deseq2_clipped = log2FoldChange_deseq2 %>%
          dplyr::case_when(
            . >= cutoff_logfc ~ cutoff_logfc,
            . <= cutoff_logfc_neg ~ cutoff_logfc_neg,
            .default = .
          ),
        logFC_edger_clipped = logFC_edger %>%
          dplyr::case_when(
            . >= cutoff_logfc ~ cutoff_logfc,
            . <= cutoff_logfc_neg ~ cutoff_logfc_neg,
            .default = .
          ),
        logFC_limma_clipped = logFC_limma %>%
          dplyr::case_when(
            . >= cutoff_logfc ~ cutoff_logfc,
            . <= cutoff_logfc_neg ~ cutoff_logfc_neg,
            .default = .
          )
      ) %>%
      # Clip padj to the threshold
      dplyr::mutate(
        padj_deseq2_clipped = padj_deseq2 %>%
          dplyr::case_when(
            . <= cutoff_padj ~ cutoff_padj,
            .default = .
          ),
        FDR_edger_clipped = FDR_edger %>%
          dplyr::case_when(
            . <= cutoff_padj ~ cutoff_padj,
            .default = .
          ),
        adj.P.Val_limma_clipped = adj.P.Val_limma %>%
          dplyr::case_when(
            . <= cutoff_padj ~ cutoff_padj,
            .default = .
          )
      ) %>%
      # Add custom shapes to dots to identify clipped genes
      ## Normal dots: Shape 19
      ## Clipped (positive): Shape -9658
      ## Clipped (negative): Shape -9668
      ## Add names for legend
      dplyr::mutate(
        volcano_shape_deseq2 = dplyr::case_when(
          .$log2FoldChange_deseq2 >= cutoff_logfc ~ -9658,
          .$log2FoldChange_deseq2 <= cutoff_logfc_neg ~ -9668,
          .$padj_deseq2 <= cutoff_padj ~ 17,
          .default = 19
        ),
        volcano_shape_edger = dplyr::case_when(
          .$logFC_edger >= cutoff_logfc ~ -9658,
          .$logFC_edger <= cutoff_logfc_neg ~ -9668,
          .$FDR_edger <= cutoff_padj ~ 17,
          .default = 19
        ),
        volcano_shape_limma = dplyr::case_when(
          .$logFC_limma >= cutoff_logfc ~ -9658,
          .$logFC_limma <= cutoff_logfc_neg ~ -9668,
          .$adj.P.Val_limma <= cutoff_padj ~ 17,
          .default = 19
        )
      ) %>%
      # Add custom size
      dplyr::mutate(
        volcano_size_deseq2 = dplyr::case_when(
          .$log2FoldChange_deseq2 >= cutoff_logfc ~ 3,
          .$log2FoldChange_deseq2 <= cutoff_logfc_neg ~ 3,
          .$padj_deseq2 <= cutoff_padj ~ 3,
          .default = 1
        ),
        volcano_size_edger = dplyr::case_when(
          .$logFC_edger >= cutoff_logfc ~ 3,
          .$logFC_edger <= cutoff_logfc_neg ~ 3,
          .$FDR_edger <= cutoff_padj ~ 3,
          .default = 1
        ),
        volcano_size_limma = dplyr::case_when(
          .$logFC_limma >= cutoff_logfc ~ 3,
          .$logFC_limma <= cutoff_logfc_neg ~ 3,
          .$adj.P.Val_limma <= cutoff_padj ~ 3,
          .default = 1
        )
      )
    # Set names for shapes
    names(results$volcano_shape_deseq2) <- dplyr::case_when(
      results$volcano_shape_deseq2 == -9658 ~ stringr::str_c(
        "logFC >",
        cutoff_logfc
      ),
      results$volcano_shape_deseq2 == -9668 ~ stringr::str_c(
        "logFC < -",
        cutoff_logfc
      ),
      results$volcano_shape_deseq2 == 17 ~ stringr::str_c(
        "padj <",
        cutoff_padj
      ),
      results$volcano_shape_deseq2 == 19 ~ "Unclipped"
    )

    names(results$volcano_shape_edger) <- dplyr::case_when(
      results$volcano_shape_edger == -9658 ~ stringr::str_c(
        "logFC >",
        cutoff_logfc
      ),
      results$volcano_shape_edger == -9668 ~ stringr::str_c(
        "logFC < -",
        cutoff_logfc
      ),
      results$volcano_shape_edger == 17 ~ stringr::str_c("padj <", cutoff_padj),
      results$volcano_shape_edger == 19 ~ "Unclipped"
    )

    names(results$volcano_shape_limma) <- dplyr::case_when(
      results$volcano_shape_limma == -9658 ~ stringr::str_c(
        "logFC >",
        cutoff_logfc
      ),
      results$volcano_shape_limma == -9668 ~ stringr::str_c(
        "logFC < -",
        cutoff_logfc
      ),
      results$volcano_shape_limma == 17 ~ stringr::str_c("padj <", cutoff_padj),
      results$volcano_shape_limma == 19 ~ "Unclipped"
    )

    # Return data
    return(results)
  })


# Create volcano plots
plots_volcano <- top_joined_volcano %>%
  purrr::imap(
    \(x, idx) {
      # Create DESeq2 volcano
      plot_deseq2 <- EnhancedVolcano::EnhancedVolcano(
        x,
        lab = x$gene_name,
        x = "log2FoldChange_deseq2_clipped",
        y = "padj_deseq2_clipped",
        xlim = c(cutoff_logfc * -1, cutoff_logfc),
        ylim = c(0, -log10(cutoff_padj)),
        xlab = bquote(~ Log[2] ~ "fold change (DESeq2)"),
        ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
        axisLabSize = 8,
        title = str_c(idx, " DESeq2"),
        titleLabSize = 8,
        subtitle = NULL,
        caption = NULL,
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = x$volcano_size_deseq2,
        labSize = 2,
        boxedLabels = FALSE,
        shapeCustom = x$volcano_shape_deseq2,
        legendPosition = "none",
        drawConnectors = TRUE,
        widthConnectors = 0.1,
        colConnectors = "grey60",
        arrowheads = FALSE,
        gridlines.major = FALSE,
        gridlines.minor = FALSE
      )

      # Create edgeR volcano
      plot_edger <- EnhancedVolcano::EnhancedVolcano(
        x,
        lab = x$gene_name,
        x = "logFC_edger_clipped",
        y = "FDR_edger_clipped",
        xlim = c(cutoff_logfc * -1, cutoff_logfc),
        ylim = c(0, -log10(cutoff_padj)),
        xlab = bquote(~ Log[2] ~ "fold change (edgeR)"),
        ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
        axisLabSize = 8,
        title = str_c(idx, " edgeR"),
        titleLabSize = 8,
        subtitle = NULL,
        caption = NULL,
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = x$volcano_size_edger,
        labSize = 2,
        boxedLabels = FALSE,
        shapeCustom = x$volcano_shape_edger,
        legendPosition = "none",
        drawConnectors = TRUE,
        widthConnectors = 0.1,
        colConnectors = "grey60",
        arrowheads = FALSE,
        gridlines.major = FALSE,
        gridlines.minor = FALSE
      )

      # Create limma volcano
      plot_limma <- EnhancedVolcano::EnhancedVolcano(
        x,
        lab = x$gene_name,
        x = "logFC_limma_clipped",
        y = "adj.P.Val_limma_clipped",
        xlim = c(cutoff_logfc * -1, cutoff_logfc),
        ylim = c(0, -log10(cutoff_padj)),
        xlab = bquote(~ Log[2] ~ "fold change (limma)"),
        ylab = bquote(~ -Log[10] ~ "adjusted p-value"),
        axisLabSize = 8,
        title = str_c(idx, " limma"),
        titleLabSize = 8,
        subtitle = NULL,
        caption = NULL,
        pCutoff = 0.05,
        FCcutoff = 1,
        pointSize = x$volcano_size_limma,
        labSize = 2,
        boxedLabels = FALSE,
        shapeCustom = x$volcano_shape_limma,
        legendPosition = "none",
        drawConnectors = TRUE,
        widthConnectors = 0.1,
        colConnectors = "grey60",
        arrowheads = FALSE,
        gridlines.major = FALSE,
        gridlines.minor = FALSE
      )

      # Return data
      plots_out <- list(
        deseq2 = plot_deseq2,
        edger = plot_edger,
        limma = plot_limma
      )
      return(plots_out)
    }
  ) %>%
  unlist(recursive = FALSE)

grid_volcano <- patchwork::wrap_plots(plots_volcano)

# Create sets of Venn diagrams
plots_venn <- top_joined %>%
  purrr::imap(\(.toptable, .contrast_name) {
    list_genes <- list(
      DESeq2 = .toptable %>%
        dplyr::filter(padj_deseq2 <= 0.05) %>%
        .$ensembl_id,
      edgeR = .toptable %>%
        dplyr::filter(FDR_edger <= 0.05) %>%
        .$ensembl_id,
      limma = .toptable %>%
        dplyr::filter(adj.P.Val_limma <= 0.05) %>%
        .$ensembl_id
    )

    plot <- list_genes %>%
      ggVennDiagram::ggVennDiagram() +
      colorspace::scale_fill_continuous_sequential(palette = "Purple-Yellow") +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = .2)) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(face = "bold")
      ) +
      ggplot2::labs(
        title = .contrast_name,
        fill = "DEGs"
      )

    # Return data
    return(plot)
  })

grid_venn <- plots_venn %>%
  patchwork::wrap_plots(nrow = 1)

# Save data
## Tables
openxlsx::write.xlsx(
  top_joined,
  file = here::here(
    "output",
    "data_expression",
    "pilot_dge.xlsx"
  ),
  asTable = TRUE
)

## Volcano plots
ggplot2::ggsave(
  here::here(
    "output",
    "plots_QC",
    "DGE pilot.png"
  ),
  grid_volcano,
  width = 16,
  height = 9,
  scale = 1
)

## Venn diagrams
ggplot2::ggsave(
  here::here(
    "output",
    "plots_QC",
    "DGE pilot venn.png"
  ),
  grid_venn,
  width = 8,
  height = 4,
  scale = 1.2
)
