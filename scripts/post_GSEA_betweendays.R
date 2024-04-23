# Run GSEA
# Run prereq scripts from a master file first
# env_prep
# dge_cellpops_as_fixed
# dge_extract_interaction_genelist
# post_GSEA_function

# GSEA between days of passage

run_GSEA(genes_D5vsD3_P5_TvsUT,
         "D5vsD3_P5_TvsUT")

run_GSEA(genes_D5vsD3_P7_TvsUT,
         "D5vsD3_P7_TvsUT")

run_GSEA(genes_D5vsD3_P13_TvsUT,
         "D5vsD3_P13_TvsUT")

# GSEA between passages, day 3

run_GSEA(genes_P7vsP5_D3_TvsUT,
         'P7vsP5_D3_TvsUT')

run_GSEA(genes_P13vsP7_D3_TvsUT,
         'P13vsP7_D3_TvsUT')

run_GSEA(genes_P13vsP5_D3_TvsUT,
         'P13vsP5_D3_TvsUT')

# GSEA between passages, day 5

run_GSEA(genes_P7vsP5_D5_TvsUT,
         'P7vsP5_D5_TvsUT')

run_GSEA(genes_P13vsP7_D5_TvsUT,
         'P13vsP7_D5_TvsUT')

run_GSEA(genes_P13vsP5_D5_TvsUT,
         'P13vsP5_D5_TvsUT')