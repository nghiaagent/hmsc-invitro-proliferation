# Run GSEA
# Run prereq scripts from a master file first
# env_prep
# dge_cellpops_as_fixed
# dge_extract_betweenpassages_genelist
# post_GSEA_function

# GSEA between passages, at each treatment and day

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

# GSEA between days, at each treatment and passage

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

# GSEA interaction term

## Effect of treatment over days

run_GSEA(genes_D5vsD3_P5_TvsUT,
         "D5vsD3_P5_TvsUT")

run_GSEA(genes_D5vsD3_P7_TvsUT,
         "D5vsD3_P7_TvsUT")

run_GSEA(genes_D5vsD3_P13_TvsUT,
         "D5vsD3_P13_TvsUT")

## Effect of treatment over passages

### At D3

run_GSEA(genes_P7vsP5_D3_TvsUT,
         'P7vsP5_D3_TvsUT')

run_GSEA(genes_P13vsP7_D3_TvsUT,
         'P13vsP7_D3_TvsUT')

run_GSEA(genes_P13vsP5_D3_TvsUT,
         'P13vsP5_D3_TvsUT')

### At D5

run_GSEA(genes_P7vsP5_D5_TvsUT,
         'P7vsP5_D5_TvsUT')

run_GSEA(genes_P13vsP7_D5_TvsUT,
         'P13vsP7_D5_TvsUT')

run_GSEA(genes_P13vsP5_D5_TvsUT,
         'P13vsP5_D5_TvsUT')