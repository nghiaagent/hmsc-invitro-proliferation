# Run GSEA
# Run prereq scripts from a master file first
# env_prep
# dge_cellpops_as_fixed
# dge_extract_interaction_genelist
# post_GSEA_function

# GSEA between days of passage

## UT

run_GSEA(genes_D5vsD3_UT_P5,
         "D5vsD3_UT_P5")

run_GSEA(genes_D5vsD3_UT_P7,
         "D5vsD3_UT_P7")

run_GSEA(genes_D5vsD3_UT_P13,
         "D5vsD3_UT_P13")

## T

run_GSEA(genes_D5vsD3_T_P5,
         "D5vsD3_T_P5")

run_GSEA(genes_D5vsD3_T_P7,
         "D5vsD3_T_P7")

run_GSEA(genes_D5vsD3_T_P13,
         "D5vsD3_T_P13")

