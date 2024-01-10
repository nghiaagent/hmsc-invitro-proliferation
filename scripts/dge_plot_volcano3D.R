# Load data

source("./scripts/dge_cellpops_as_fixed.R")

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
  topTable(fit_contrasts,           number = Inf, sort.by = "none")$P.Value,
  topTable(fit_contrasts, coef = 4, number = Inf, sort.by = "none")$P.Value,
  topTable(fit_contrasts, coef = 5, number = Inf, sort.by = "none")$P.Value,
  topTable(fit_contrasts, coef = 6, number = Inf, sort.by = "none")$P.Value
)

polar_padj <- cbind(
  topTable(fit_contrasts,           number = Inf, sort.by = "none")$adj.P.Val,
  topTable(fit_contrasts, coef = 4, number = Inf, sort.by = "none")$adj.P.Val,
  topTable(fit_contrasts, coef = 5, number = Inf, sort.by = "none")$adj.P.Val,
  topTable(fit_contrasts, coef = 6, number = Inf, sort.by = "none")$adj.P.Val
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
              r_axis_ticks = c(0.5,
                               1,
                               1.5,
                              2,
                               2.5,
                               3.0
                               ),
              z_axis_ticks = c(0.5),
              axis_angle = 12/6)

