# Perform camera on logCPM produced by limma-voom

# List of gene sets to be used for camera
# GO (CC, BP, MF)
# KEGG
# Reactome
# MSigDB h, c2/CGP, c3

#### Define function for converting GMT files to limma::camera format ####

get_camera_gs <- function(con) {
  getGmt(con = con) %>%
    geneIds() %>%
    ids2indices(identifiers = fit_contrasts$genes$ENTREZID)
}

#### Obtain gene sets ####

# GO

## GOBP
## Represented by MSigDB c5/GO/BP

msigdb_GOBP <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.entrez.gmt")

## GOMF
## Represented by MSigDB c5/GO/MF

msigdb_GOMF <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.mf.v2023.2.Hs.entrez.gmt")

## GOCC
## Represented by MSigDB c5/GO/CC

msigdb_GOCC <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c5.go.cc.v2023.2.Hs.entrez.gmt")

# MSigDB

msigdb_h  <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.entrez.gmt")

msigdb_c2 <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cgp.v2023.2.Hs.entrez.gmt")

msigdb_c3 <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c3.all.v2023.2.Hs.entrez.gmt")

# ReactomePA
# Represented by MSigDB c2/CP/Reactome

msigdb_reactome <-
  get_camera_gs("./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cp.reactome.v2023.2.Hs.entrez.gmt")

# KEGG
# Represented by MSigDB c2/CP/KEGG_LEGACY

msigdb_KEGG <-
  get_camera_gs(
    "./input/genesets/msigdb_v2023.2.Hs_GMTs/c2.cp.kegg_legacy.v2023.2.Hs.entrez.gmt"
  )

#### Define camera helper function ####
# Function runs camera on MArrayLM fit object provided by voomLmFit,
# allowing selection of contrast coefficient and threshold
# Export files to relevant folder

run_camera <- function(fit = NULL, coef = NULL, inter.gene.cor = 0.01, sort = TRUE) {
  
  # Check for errors in input
  ## Check errors in fit object
  
  if (!class(fit) == "MArrayLM" | is.null(fit$EList) | is.null(fit$contrasts)) {
    stop("Fit object is not a contrasted MArrayLM object provided by voomLmFit.\nIs keep.EList = TRUE? Have you run contrasts.fit after voomLmFit?")
  }
  
  ## Check errors in inter-gene correlation
  
  if (is.null(inter.gene.cor) == TRUE) {
    
    message(str_c("inter-gene correlation estimated by camera"))
    
  }
  else {
    
    if ((is.numeric(inter.gene.cor) == FALSE |
         abs(inter.gene.cor) > 1)) {
      stop("inter-gene correlation must be between -1 and 1")
    }
    message(str_c("inter-gene correlation set at ", inter.gene.cor))
    
  }
  
  ## Check errors in sort option
  
  if (class(sort) != "logical") {
    stop("sort selection must be TRUE or FALSE")
  }
  
  # Set up output folder
  
  name_output <- as.character(colnames(fit$coefficients)[coef])
  
  path_output <- file.path(getwd(), 
                           'output',
                           'data_enrichment', 
                           'camera', 
                           name_output)
  
  message(str_c("Output camera results to", path_output, 
                sep = " "))
  
  
  if (!dir.exists(path_output)) {
    
    dir.create(path_output, recursive = TRUE)
    
  }
  
  # Get design matrix
  
  design <- fit$design
  
  # Run camera on GO gene sets
  
  message(str_c("Running GO enrichment for", name_output, sep = " "))
  
  camera_GOBP <- camera(
    fit$EList$E,
    msigdb_GOBP,
    fit$design,
    fit$contrasts[, coef],
    sort = sort,
    inter.gene.cor = inter.gene.cor
  ) %>%
    mutate(name = rownames(.)) 
  
  camera_GOCC <- camera(
    fit$EList$E,
    msigdb_GOCC,
    fit$design,
    fit$contrasts[, coef],
    sort = sort,
    inter.gene.cor = inter.gene.cor
  ) %>%
    mutate(name = rownames(.)) 
  
  camera_GOMF <- camera(
    fit$EList$E,
    msigdb_GOMF,
    fit$design,
    fit$contrasts[, coef],
    sort = sort,
    inter.gene.cor = inter.gene.cor
  ) %>%
    mutate(name = rownames(.)) 
  
  message(str_c("Running MSigDB enrichment for", name_output, sep = " "))
  
  camera_h <- camera(
    fit$EList$E,
    msigdb_h,
    fit$design,
    fit$contrasts[, coef],
    sort = sort,
    inter.gene.cor = inter.gene.cor
  ) %>%
    mutate(name = rownames(.)) 
  
  camera_c2 <- camera(
    fit$EList$E,
    msigdb_c2,
    fit$design,
    fit$contrasts[, coef],
    sort = sort,
    inter.gene.cor = inter.gene.cor
  ) %>%
    mutate(name = rownames(.)) 
  
  camera_c3 <- camera(
    fit$EList$E,
    msigdb_c3,
    fit$design,
    fit$contrasts[, coef],
    sort = sort,
    inter.gene.cor = inter.gene.cor
  ) %>%
    mutate(name = rownames(.)) 
  
  message(str_c("Running Reactome enrichment for", name_output, sep = " "))
  
  camera_reactome <- camera(
    fit$EList$E,
    msigdb_reactome,
    fit$design,
    fit$contrasts[, coef],
    sort = sort,
    inter.gene.cor = inter.gene.cor
  ) %>%
    mutate(name = rownames(.)) 
  
  message(str_c("Running KEGG enrichment for", name_output, sep = " "))
  
  camera_KEGG <- camera(
    fit$EList$E,
    msigdb_KEGG,
    fit$design,
    fit$contrasts[, coef],
    sort = sort,
    inter.gene.cor = inter.gene.cor
  ) %>%
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
  writexl::write_xlsx(camera_list, 
                      path = file.path(path_output, "camera_results.xlsx"))
  return(camera_list)
  
}
