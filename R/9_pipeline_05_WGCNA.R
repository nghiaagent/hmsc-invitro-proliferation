here::i_am("R/9_pipeline_05_WGCNA.R")

########################
# Run WGCNA
########################

# Run scripts of pipeline
c(
  "5_post_WGCNA_1_load_and_QC.R",
  "5_post_WGCNA_2_construct_net.R",
  "5_post_WGCNA_3_find_hubs.R",
  "5_post_WGCNA_4_extract_module_genes.R",
  "5_post_WGCNA_5_export_Cytoscape.R"
) %>%
  walk(
    \(x) {
      message(paste0(
        "Sourcing ",
        here::here("scripts", x)
      ))
      source(here::here("R", x), echo = TRUE, verbose = FALSE)
    },
    .progress = TRUE
  )
