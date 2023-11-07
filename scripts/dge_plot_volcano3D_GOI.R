# Load data

source("./scripts/dge_selected_genes.R")

# Supply p-values to the 3D volcano plot
## First column: ANOVA p-values
## Remaining columns: Comparisons
## 2: A vs B (P5 vs P7)
## 3: A vs C (P5 vs P13)
## 4: B vs C (P7 vs P13)

table_design$condition_ID <- factor(table_design$condition_ID,
          levels = c("P5D3Untreated",
                     "P7D3Untreated",
                     "P13D3Untreated"),
          labels = c("P5",
                     "P7",
                     "P13"))

polar_pvals <- cbind(
  topTable(fit_small,           number = Inf, sort.by = "none")$P.Value,
  topTable(fit_small, coef = 1, number = Inf, sort.by = "none")$P.Value,
  topTable(fit_small, coef = 2, number = Inf, sort.by = "none")$P.Value,
  topTable(fit_small, coef = 3, number = Inf, sort.by = "none")$P.Value
)

polar_padj <- cbind(
  topTable(fit_small,           number = Inf, sort.by = "none")$adj.P.Val,
  topTable(fit_small, coef = 1, number = Inf, sort.by = "none")$adj.P.Val,
  topTable(fit_small, coef = 2, number = Inf, sort.by = "none")$adj.P.Val,
  topTable(fit_small, coef = 3, number = Inf, sort.by = "none")$adj.P.Val
)

## Construct volcano3d object

quant_DGE_voom_small <- quant_DGE_voom[quant_DGE_voom$genes$ENTREZID %in% entrezid_small,]

rownames(quant_DGE_voom_small$E) <- quant_DGE_voom_small$genes$GENENAME

polar_manual <- polar_coords(
  outcome = table_design$condition_ID,
  data = t(quant_DGE_voom_small$E),
  pvals = polar_pvals,
  padj = polar_padj
)

rownames(polar_manual@pvals) <- quant_DGE_voom_small$genes$GENENAME
colnames(polar_manual@pvals) <- c("ANOVA", "P7vsP5", "P13vsP5", "P13vsP7")
rownames(polar_manual@padj) <- quant_DGE_voom_small$genes$GENENAME
colnames(polar_manual@padj) <- c("ANOVA", "P7vsP5", "P13vsP5", "P13vsP7")

## Plot

radial_plotly(polar_manual,
              axis_angle = 0.8,
              label_rows = c("GPC1",
                             "GPC4",
                             "SDC1",
                             "SDC4"
                             ),
              label_size = 28,
              arrow_length = 120)

## 3D Volcano plot
## Trio boxplot of GOIs

boxplot_trio(polar = polar_manual,
             value = "NES")
