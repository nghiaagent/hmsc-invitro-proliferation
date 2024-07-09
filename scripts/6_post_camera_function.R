# Perform camera on t-statistics produced by limma-voom
# Load data

quant_DGE_voom <- readRDS("./output/data_expression/post_DGE/quant_DGE_voom.RDS")
fit_contrasts <- readRDS("./output/data_expression/post_DGE/fit_contrasts.RDS")

# List of gene sets to be used for camera
# GO (CC, BP, MF)
# KEGG
# Reactome
# MSigDB h, c2/CGP, c3

# Define function for converting GMT files to limma::camera gene sets format

get_camera_gs <- function(con) {
  getGmt(con = con) %>%
    geneIds() %>%
    ids2indices(identifiers = fit_contrasts$genes$ENTREZID)
}

# Obtain gene sets

## GO

### GOBP
### Represented by MSigDB c5/GO/BP

msigdb_GOBP <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.entrez.gmt")

### GOMF
### Represented by MSigDB c5/GO/MF

msigdb_GOMF <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.mf.v2023.2.Hs.entrez.gmt")

### GOCC
### Represented by MSigDB c5/GO/CC

msigdb_GOCC <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.cc.v2023.2.Hs.entrez.gmt")

## MSigDB

msigdb_h  <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.entrez.gmt")

msigdb_c2 <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cgp.v2023.2.Hs.entrez.gmt")

msigdb_c3 <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c3.all.v2023.2.Hs.entrez.gmt")

## ReactomePA
## Represented by MSigDB c2/CP/Reactome

msigdb_reactome <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cp.reactome.v2023.2.Hs.entrez.gmt")

## KEGG
## Represented by MSigDB c2/CP/KEGG_LEGACY

msigdb_KEGG <-
  get_camera_gs(
    "./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cp.kegg_legacy.v2023.2.Hs.entrez.gmt"
  )

# Test camera

## CameraPR on t statistics

run_camera <- function(fit, coefficient) {
  name_output <- as.character(colnames(matrix_contrasts)[coefficient])
  
  path_output <- file.path(getwd(), 'output', 'data_enrichment', 'camera', name_output)
  
  message(str_c("Output camera results to", path_output, sep = " "))
  
  if (!dir.exists(path_output)) {
    dir.create(path_output, recursive = TRUE)
  }
  
  # Run camera on GO gene sets
  
  message(str_c("Running GO enrichment for", name_output, sep = " "))
  
  camera_GOBP <- camera(quant_DGE_voom$E,
                        msigdb_GOBP,
                        design,
                        matrix_contrasts[, coefficient],
                        sort = FALSE,
                        inter.gene.cor = NULL) %>%
    mutate(name = rownames(.))
  
  camera_GOCC <- camera(quant_DGE_voom$E,
                        msigdb_GOCC,
                        design,
                        matrix_contrasts[, coefficient],
                        sort = FALSE,
                        inter.gene.cor = NULL) %>%
    mutate(name = rownames(.))
  
  camera_GOMF <- camera(quant_DGE_voom$E,
                        msigdb_GOMF,
                        design,
                        matrix_contrasts[, coefficient],
                        sort = FALSE,
                        inter.gene.cor = NULL) %>%
    mutate(name = rownames(.))
  
  message(str_c("Running MSigDB enrichment for", name_output, sep = " "))
  
  camera_h <- camera(quant_DGE_voom$E,
                     msigdb_h,
                     design,
                     matrix_contrasts[, coefficient],
                     sort = FALSE,
                     inter.gene.cor = NULL) %>%
    mutate(name = rownames(.))
  
  camera_c2 <- camera(quant_DGE_voom$E,
                      msigdb_c2,
                      design,
                      matrix_contrasts[, coefficient],
                      sort = FALSE,
                      inter.gene.cor = NULL) %>%
    mutate(name = rownames(.))
  
  camera_c3 <- camera(quant_DGE_voom$E,
                      msigdb_c3,
                      design,
                      matrix_contrasts[, coefficient],
                      sort = FALSE,
                      inter.gene.cor = NULL) %>%
    mutate(name = rownames(.))
  
  message(str_c("Running Reactome enrichment for", name_output, sep = " "))
  
  camera_reactome <- camera(quant_DGE_voom$E,
                            msigdb_reactome,
                            design,
                            matrix_contrasts[, coefficient],
                            sort = FALSE,
                            inter.gene.cor = NULL) %>%
    mutate(name = rownames(.))
  
  message(str_c("Running KEGG enrichment for", name_output, sep = " "))
  
  camera_KEGG <- camera(quant_DGE_voom$E,
                        msigdb_KEGG,
                        design, 
                        matrix_contrasts[, coefficient],
                        sort = FALSE,
                        inter.gene.cor = NULL) %>%
    mutate(name = rownames(.))
  
  camera_list <- list(
    "GOBP" = camera_GOBP,
    "GOCC" = camera_GOCC,
    "GOMF" = camera_GOMF,
    "KEGG" = camera_KEGG,
    "MSigDB_h" = camera_h,
    "MSigDB_c2" = camera_c2,
    "MSigDB_c3" = camera_c3,
    "Reactome" = camera_reactome
  )
  
  saveRDS(camera_list, file = file.path(path_output, "camera_results.RDS"))
  
  return(camera_list)
  
}
