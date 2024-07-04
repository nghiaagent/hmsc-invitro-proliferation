# Prepare dataset.

fit_contrasts <-
  readRDS(file = "./output/data_expression/post_DGE/fit_contrasts.RDS")

list_top <- map(1:ncol(fit_contrasts$coefficients),
         \ (x) topTable(fit_contrasts, coef = x, n = Inf))

names(list_top) <- colnames(fit_contrasts$contrasts)

writexl::write_xlsx(list_top, path = file.path("output",
                                               "data_expression",
                                               "post_DGE",
                                               "genes_all.xlsx"))