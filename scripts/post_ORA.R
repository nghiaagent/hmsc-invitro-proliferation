# Run ORA
# Run prereq scripts from a master file first
# env_prep
# dge_cellpops_as_fixed
# dge_extract_venn_genes
# post_ORA_function

# ORA between days of passage

run_ORA(genes_D3_Untreated_P13vsP5_and_P13vsP7_only,
        "D3_Untreated_P13vsP5_and_P13vsP7_only")
# 
# run_ORA(genes_D3_Untreated_P13vsP5_only,
#         "D3_Untreated_P13vsP5_only")

run_ORA(genes_D3_Untreated_P13vsP7_only,
        "D3_Untreated_P13vsP7_only")

run_ORA(genes_D3_Untreated_P7vsP5_only,
        "D3_Untreated_P7vsP5_only")

run_ORA(genes_D5_Untreated_P13vsP5_and_P13vsP7_only,
        "D5_Untreated_P13vsP5_and_P13vsP7_only")

run_ORA(genes_D5_Untreated_P13vsP5_only,
        "D5_Untreated_P13vsP5_only")

run_ORA(genes_D5_Untreated_P13vsP7_only,
        "D5_Untreated_P13vsP7_only")

run_ORA(genes_D5_Untreated_P7vsP5_only,
        "D5_Untreated_P7vsP5_only")

# ORA days x treat

run_ORA(genes_D5vsD3_P13_Overlap,
        "D5vsD3_P13_Overlap")

run_ORA(genes_D5vsD3_P13_Untreated_only,
        "D5vsD3_P13_Untreated_only")

run_ORA(genes_D5vsD3_P5_Overlap,
        "D5vsD3_P5_Overlap")

run_ORA(genes_D5vsD3_P5_Treated_only,
        "D5vsD3_P5_Treated_only")

run_ORA(genes_D5vsD3_P5_Untreated_only,
        "D5vsD3_P5_Untreated_only")

run_ORA(genes_D5vsD3_P7_Overlap,
        "D5vsD3_P7_Overlap")

run_ORA(genes_D5vsD3_P7_Treated_only,
        "D5vsD3_P7_Treated_only")

run_ORA(genes_D5vsD3_P7_Untreated_only,
        "D5vsD3_P7_Untreated_only")

# ORA passage x treat at each day

## D3

run_ORA(genes_P13vsP5_D3_Overlap,
        "P13vsP5_D3_Overlap")

run_ORA(genes_P13vsP7_D3_Overlap,
        "P13vsP7_D3_Overlap")

run_ORA(genes_P7vsP5_D3_Overlap,
        "P7vsP5_D3_Overlap")

run_ORA(genes_P13vsP5_D3_Treated_only,
        "P13vsP5_D3_Treated_only")

run_ORA(genes_P13vsP7_D3_Treated_only,
        "P13vsP7_D3_Treated_only")

run_ORA(genes_P7vsP5_D3_Treated_only,
        "P7vsP5_D3_Treated_only")

run_ORA(genes_P13vsP5_D3_Untreated_only,
        "P13vsP5_D3_Untreated_only")

run_ORA(genes_P13vsP7_D3_Untreated_only,
        "P13vsP7_D3_Untreated_only")

## D5

run_ORA(genes_P13vsP5_D5_Overlap,
        "P13vsP5_D5_Overlap")

run_ORA(genes_P13vsP7_D5_Overlap,
        "P13vsP7_D5_Overlap")

run_ORA(genes_P7vsP5_D5_Overlap,
        "P7vsP5_D5_Overlap")

run_ORA(genes_P13vsP5_D5_Untreated_only,
        "P13vsP5_D5_Untreated_only")

run_ORA(genes_P13vsP7_D5_Untreated_only,
        "P13vsP7_D5_Untreated_only")

run_ORA(genes_P7vsP5_D5_Untreated_only,
        "P7vsP5_D5_Untreated_only")