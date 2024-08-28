# Run scripts of limma::voom pipeline

list_files <- c(
  "0_meta_env_prep.R",
  "1_ensembl_voom_pre_quant_import.R",
  "1_ensembl_voom_pre_quality_control_filter.R",
  "2_ensembl_voom_dge.R"
)

walk(list_files, \ (x) source(here::here("scripts", x)))