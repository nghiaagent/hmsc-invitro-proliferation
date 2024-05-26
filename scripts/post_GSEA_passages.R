# Run GSEA
# Run prereq scripts from a master file first
# env_prep
# dge_cellpops_as_fixed
# dge_extract_betweenpassages_genelist
# post_GSEA_function

# GSEA between passaged

## UT D3

run_GSEA(genes_P7vsP5_UT_D3,
         "P7vsP5_UT_D3")

run_GSEA(genes_P13vsP7_UT_D3,
         "P13vsP7_UT_D3")

run_GSEA(genes_P13vsP5_UT_D3,
         "P13vsP5_UT_D3")

## T D3

run_GSEA(genes_P7vsP5_T_D3,
         "P7vsP5_T_D3")

run_GSEA(genes_P13vsP7_T_D3,
         "P13vsP7_T_D3")

run_GSEA(genes_P13vsP5_T_D3,
         "P13vsP5_T_D3")

## UT D5

run_GSEA(genes_P7vsP5_UT_D5,
         "P7vsP5_UT_D5")

run_GSEA(genes_P13vsP7_UT_D5,
         "P13vsP7_UT_D5")

run_GSEA(genes_P13vsP5_UT_D5,
         "P13vsP5_UT_D5")

## T D5

run_GSEA(genes_P7vsP5_T_D5,
         "P7vsP5_T_D5")

run_GSEA(genes_P13vsP7_T_D5,
         "P13vsP7_T_D5")

run_GSEA(genes_P13vsP5_T_D5,
         "P13vsP5_T_D5")
