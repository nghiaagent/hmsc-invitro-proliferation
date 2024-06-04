# Run ORA
# Run prereq scripts from a master file first
# env_prep
# dge_cellpops_as_fixed
# dge_extract_venn_genes
# post_WGCNA_extract_module_genes
# post_ORA_function

# ORA between passages

run_ORA(genes_P7vsP5_D3_Untreated$ENTREZID,
        "P7vsP5_D3_Untreated")

run_ORA(genes_P13vsP7_D3_Untreated$ENTREZID,
        "P13vsP7_D3_Untreated")

run_ORA(genes_P13vsP5_D3_Untreated$ENTREZID,
        "P13vsP5_D3_Untreated")

# ORA on overlap in changes between passages in untreated cells

run_ORA(genes_D3_Untreated_P13vsP5_and_P13vsP7_only$ENTREZID,
        "D3_Untreated_P13vsP5_and_P13vsP7_only")

run_ORA(genes_D3_Untreated_P13vsP5_only$ENTREZID,
        "D3_Untreated_P13vsP5_only")

run_ORA(genes_D3_Untreated_P13vsP7_only$ENTREZID,
        "D3_Untreated_P13vsP7_only")

run_ORA(genes_D3_Untreated_P7vsP5_only$ENTREZID,
        "D3_Untreated_P7vsP5_only")

run_ORA(genes_D5_Untreated_P13vsP5_and_P13vsP7_only$ENTREZID,
        "D5_Untreated_P13vsP5_and_P13vsP7_only")

run_ORA(genes_D5_Untreated_P13vsP5_only$ENTREZID,
        "D5_Untreated_P13vsP5_only")

run_ORA(genes_D5_Untreated_P13vsP7_only$ENTREZID,
        "D5_Untreated_P13vsP7_only")

run_ORA(genes_D5_Untreated_P7vsP5_only$ENTREZID,
        "D5_Untreated_P7vsP5_only")

# ORA on GCN gene sets
## Modules correlated with passage / day selected
## And modules differentially regulated with GSVA

## WGCNA with all samples - Modules: Red, Turquoise

### Corr with trait - Modules: Red, Turquoise

run_ORA(GCN_genelists$WGCNA_allsamples_turquoise,
        "WGCNA_allsamples_turquoise")

run_ORA(GCN_genelists$WGCNA_allsamples_red,
        "WGCNA_allsamples_red")

### DEM in GSVA - Modules: Black, Green, Magenta

run_ORA(GCN_genelists$WGCNA_allsamples_black,
        "WGCNA_allsamples_black")

run_ORA(GCN_genelists$WGCNA_allsamples_green,
        "WGCNA_allsamples_green")

run_ORA(GCN_genelists$WGCNA_allsamples_magenta,
        "WGCNA_allsamples_magenta")

## WGCNA with D3 UT samples 
### Corr with trait - Modules: Turquoise, Green, Yellow, Brown

run_ORA(GCN_genelists$WGCNA_subset_turquoise,
        "WGCNA_subset_turquoise")

run_ORA(GCN_genelists$WGCNA_subset_green,
        "WGCNA_subset_green")

run_ORA(GCN_genelists$WGCNA_subset_yellow,
        "WGCNA_subset_yellow")

run_ORA(GCN_genelists$WGCNA_subset_brown,
        "WGCNA_subset_brown")

### DEM in GSVA - Modules: Red

run_ORA(GCN_genelists$WGCNA_subset_red,
        "WGCNA_subset_red")

## CEMiTool: All 5

run_ORA(GCN_genelists$CEMiTool_M1,
        "CEMiTool_M1")

run_ORA(GCN_genelists$CEMiTool_M2,
        "CEMiTool_M2")

run_ORA(GCN_genelists$CEMiTool_M3,
        "CEMiTool_M3")

run_ORA(GCN_genelists$CEMiTool_M4,
        "CEMiTool_M4")

run_ORA(GCN_genelists$CEMiTool_M5,
        "CEMiTool_M5")