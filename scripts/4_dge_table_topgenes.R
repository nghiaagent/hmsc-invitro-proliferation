# Extract logFC and p-vals for comparisons between passages at D3 and D5

top_P7vsP5_D3 <- topTable(fit_contrasts, coef = 4,
                          number = Inf,
                          sort.by = "p")

top_P7vsP5_D5 <- topTable(fit_contrasts, coef = 13,
                          number = Inf,
                          sort.by = "p")

top_P13vsP5_D3 <- topTable(fit_contrasts, coef = 5,
                           number = Inf,
                           sort.by = "p")

top_P13vsP5_D5 <- topTable(fit_contrasts, coef = 14,
                           number = Inf,
                           sort.by = "p")

top_P13vsP7_D3 <- topTable(fit_contrasts, coef = 6,
                           number = Inf,
                           sort.by = "p")

top_P13vsP7_D5 <- topTable(fit_contrasts, coef = 15,
                           number = Inf,
                           sort.by = "p")

# Extract logFC and p-vals for Trt vs Ctrl

top_treat_P5_D3 <- topTable(fit_contrasts, coef = 1,
                          number = Inf,
                          sort.by = "p")

top_treat_P5_D5 <- topTable(fit_contrasts, coef = 10,
                          number = Inf,
                          sort.by = "p")

top_treat_P7_D3 <- topTable(fit_contrasts, coef = 2,
                           number = Inf,
                           sort.by = "p")

top_treat_P7_D5 <- topTable(fit_contrasts, coef = 11,
                           number = Inf,
                           sort.by = "p")

top_treat_P13_D3 <- topTable(fit_contrasts, coef = 3,
                           number = Inf,
                           sort.by = "p")

top_treat_P13_D5 <- topTable(fit_contrasts, coef = 12,
                           number = Inf,
                           sort.by = "p")

list_sheets <- list(
  "Phase B - A, D3" = top_P7vsP5_D3,
  "Phase C - A, D3" = top_P13vsP5_D3,
  "Phase C - B, D3" = top_P13vsP7_D3,
  "Phase B - A, D5" = top_P7vsP5_D5,
  "Phase C - A, D5" = top_P13vsP5_D5,
  "Phase C - B, D5" = top_P13vsP7_D5,
  "Treatment @ Phase A, D3" = top_treat_P5_D3,
  "Treatment @ Phase B, D3" = top_treat_P7_D3,
  "Treatment @ Phase C, D3" = top_treat_P13_D3,
  "Treatment @ Phase A, D5" = top_treat_P5_D5,
  "Treatment @ Phase B, D5" = top_treat_P7_D5,
  "Treatment @ Phase C, D5" = top_treat_P13_D5
)

write_xlsx(list_sheets, "./output/top_DEGs.xlsx")
