# Run env_prep and dge_cellpops_as_fixed in succession to load packages and run DGE

source(file.path(".",
                 "scripts",
                 "0_meta_env_prep.R"))

source(file.path(".",
                 "scripts",
                 "2_dge.R"))