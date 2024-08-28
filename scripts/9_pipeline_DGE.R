# Run scripts of limma::voom pipeline

list_files_limma <- c(
  "0_meta_env_prep.R",
  "1_ensembl_voom_pre_quant_import.R",
  "1_ensembl_voom_pre_quality_control_filter.R",
  "2_ensembl_voom_dge.R"
)


list_files_deseq2 <- c(
  "0_meta_env_prep.R",
  "1_ensembl_deseq2_pre_quant_import.R",
  "2_ensembl_deseq2_dge.R"
)

walk(list_files_limma,
     \ (x) source(here::here("scripts", x)))

walk(list_files_deseq2,
     \ (x) source(here::here("scripts", x)))