# Load data

source("./scripts/dge_cellpops_as_fixed.R")

# Supply p-values to the 3D volcano plot
## First column: ANOVA p-values
## Remaining columns: Comparisons
## 2: A vs B (P5 vs P7)
## 3: A vs C (P5 vs P13)
## 4: B vs C (P7 vs P13)

table_design$condition_ID <- factor(table_design$condition_ID,
          levels = c("P5D5Untreated",
                     "P7D5Untreated",
                     "P13D5Untreated"),
          labels = c("A",
                     "B",
                     "C"))

polar_pvals <- cbind(
  topTable(fit_contrasts,           number = Inf, sort.by = "none")$P.Value,
  topTable(fit_contrasts, coef = 13, number = Inf, sort.by = "none")$P.Value,
  topTable(fit_contrasts, coef = 14, number = Inf, sort.by = "none")$P.Value,
  topTable(fit_contrasts, coef = 15, number = Inf, sort.by = "none")$P.Value
)

polar_padj <- cbind(
  topTable(fit_contrasts,           number = Inf, sort.by = "none")$adj.P.Val,
  topTable(fit_contrasts, coef = 13, number = Inf, sort.by = "none")$adj.P.Val,
  topTable(fit_contrasts, coef = 14, number = Inf, sort.by = "none")$adj.P.Val,
  topTable(fit_contrasts, coef = 15, number = Inf, sort.by = "none")$adj.P.Val
)

## Construct volcano3d object

rownames(quant_DGE_voom$E) <- quant_DGE_voom$genes$GENENAME

polar_manual <- polar_coords(
  outcome = table_design$condition_ID,
  data = t(quant_DGE_voom$E),
  pvals = polar_pvals,
  padj = polar_padj
)

## Plot

volcano3D(polar_manual,
              axis_angle = 3/6,
              label_size = 20,
              axis_label_size = 18,
          axis_title_size = 28)

