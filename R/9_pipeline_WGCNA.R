# Run scripts of WGCNA pipeline

list_files <- c(
    "0_meta_env_prep.R",
    "1_ensembl_voom_pre_quant_import.R",
    "1_ensembl_voom_pre_quality_control_filter.R",
    "2_ensembl_voom_dge.R",
    "5_post_WGCNA_1a_load_and_QC_allsamples.R",
    "5_post_WGCNA_2a_construct_net_allsamples.R",
    "5_post_WGCNA_1b_load_and_QC_D3_UT_only.R",
    "5_post_WGCNA_2b_construct_net_D3_UT_only.R"
)

walk(list_files, \(x) source(here::here("scripts", x)))
