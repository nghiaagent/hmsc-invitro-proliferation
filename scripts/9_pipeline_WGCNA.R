# Load packages, perform DGE

source(file.path(".",
                 "scripts",
                 "env_prep.R"))

source(file.path(".",
                 "scripts",
                 "dge_cellpops_as_fixed.R"))

# Perform WGCNA with all samples

source(file.path(".",
                 "scripts",
                 "post_WGCNA_1a_load_and_QC_allsamples.R"))

source(file.path(".",
                 "scripts",
                 "post_WGCNA_2a_construct_net_allsamples.R"))

# Perform WGCNA on D3 UT samples only

source(file.path(".",
                 "scripts",
                 "post_WGCNA_1b_load_and_QC_D3_UT_only.R"))

source(file.path(".",
                 "scripts",
                 "post_WGCNA_2b_construct_net_D3_UT_only.R"))