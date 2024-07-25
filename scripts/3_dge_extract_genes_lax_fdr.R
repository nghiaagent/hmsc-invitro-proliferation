# Load data

fit_contrasts <- readRDS(file = "./output/data_expression/post_DGE/fit_contrasts.RDS")

# Extract genes differentially regulated with a relaxed FDR
# From 5% to 7.5%, 10%, 12.5%, 15%

# Calculate no. of significant genes at each FDR

## Define p-vals and contrasts to be tested

list_contrasts <- c(1:length(colnames(fit_contrasts$contrasts))) %>%
  magrittr::set_names(colnames(fit_contrasts$contrasts))

list_pvals <- c(
  test_0.05  = 0.05,
  test_0.075 = 0.075,
  test_0.1   = 0.1,
  test_0.125 = 0.125,
  test_0.15  = 0.15,
  test_0.175 = 0.175,
  test_0.2   = 0.2
)

## Get TestMatrix lists

tests_matrix <- map(.x = list_pvals,
                    .f = \(x) decideTests(fit_contrasts, p.value = x))

## Merge lists into one summary table

tests_summary <- map(.x = tests_matrix, .f = \(x) summary(x)) %>%
  do.call("rbind", .) %>%
  t() %>%
  as.data.frame() %>%
  magrittr::set_colnames(
    map(.x = list_pvals,
        .f = \(x) paste0(x, c("_down", "_ns", "_up"))) %>%
    unlist()
    )

tests_summary_compact <- tests_summary %>%
  select(map(.x = list_pvals,
             .f = \(x) paste0(x, c("_down", "_up"))) %>%
           unlist() %>%
           unname())

# Get names of significant genes at each comparison, each FDR
# Done using nested map:
# Inner map: Get names of significant genes for every contrasts
# Outer map: Performs inner map on every FDR chosen

signif_genes <- map(.x = tests_matrix,
                    .f = \(tests) {
  map(.x = list_contrasts,
      .f = \(coef) {fit_contrasts$genes[which(tests[, coef] != 0), ]$GENENAME})
}
)

# Export data

## Export summary table

openxlsx::write.xlsx(tests_summary,
                     file = file.path("output",
                                      "data_expression",
                                      "post_DGE",
                                      "FDR_screening.xlsx"),
                     asTable = TRUE)

