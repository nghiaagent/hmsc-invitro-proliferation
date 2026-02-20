# Run script to load and install packages

source("scripts/0_meta_env_prep.R", echo = TRUE, verbose = TRUE)

# Run scripts of DESeq2 pipeline

walk(
  c(
    "1_pre_quant_import.R",
    "2_dge.R",
    "4_dge_plot_volcano2D.R"
  ),
  \(x) source(here::here("scripts", x), echo = TRUE, verbose = TRUE)
)
