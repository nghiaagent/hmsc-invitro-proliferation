here::i_am("R/6_post_camera_01_prepare.R")

########################
# Prepare data for CAMERA
# Define design and contrast matrices
# Defun function to run camera
########################

# Import packages
library(DESeq2)
library(GSEABase)
library(here)
library(limma)
library(SummarizedExperiment)
library(tidyverse)

# Load data
rlog_deseq2 <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "rlog_deseq2.RDS"
  )
)

quant_deseq2 <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

# Convert rlog object to MArrayLM to pass to camera
# Method accepted by camera authors https://support.bioconductor.org/p/71534/
rlog_camera <- new(
  "EList",
  list(
    E = assay(rlog_deseq2),
    genes = rowRanges(rlog_deseq2),
    targets = colData(rlog_deseq2)
  )
)

# Define design and contrast matrices for camera
## Define design matrix
design <- model.matrix(
  design(quant_deseq2),
  data = rlog_camera$targets
)

colnames(design) <- make.names(colnames(design))

## Define contrast matrix
contrasts <- makeContrasts(
  ## Coefs 1 - 6: Treatment at each timepoint
  Trt_P5_D3 = condition_IDP5D3Treated - 0,
  Trt_P5_D5 = condition_IDP5D5Treated - condition_IDP5D5Untreated,
  Trt_P7_D3 = condition_IDP7D3Treated - condition_IDP7D3Untreated,
  Trt_P7_D5 = condition_IDP7D5Treated - condition_IDP7D5Untreated,
  Trt_P13_D3 = condition_IDP13D3Treated - condition_IDP13D3Untreated,
  Trt_P13_D5 = condition_IDP13D5Treated - condition_IDP13D5Untreated,

  ## Coefs 7 - 12: Day at each timepoint x treatment
  D5vsD3_UT_P5 = condition_IDP5D5Untreated - 0,
  D5vsD3_UT_P7 = condition_IDP7D5Untreated - condition_IDP7D3Untreated,
  D5vsD3_UT_P13 = condition_IDP13D5Untreated - condition_IDP13D3Untreated,
  D5vsD3_T_P5 = condition_IDP5D5Treated - condition_IDP5D3Treated,
  D5vsD3_T_P7 = condition_IDP7D5Treated - condition_IDP7D3Treated,
  D5vsD3_T_P13 = condition_IDP13D5Treated - condition_IDP13D3Treated,

  ## Coefs 13 - 24: Passage at each day x treatment
  P7vsP5_UT_D3 = condition_IDP7D3Untreated - 0,
  P13vsP7_UT_D3 = condition_IDP13D3Untreated - condition_IDP7D3Untreated,
  P13vsP5_UT_D3 = condition_IDP13D3Untreated - 0,
  P7vsP5_T_D3 = condition_IDP7D3Treated - condition_IDP5D3Treated,
  P13vsP7_T_D3 = condition_IDP13D3Treated - condition_IDP7D3Treated,
  P13vsP5_T_D3 = condition_IDP13D3Treated - condition_IDP5D3Treated,
  P7vsP5_UT_D5 = condition_IDP7D5Untreated - condition_IDP5D5Untreated,
  P13vsP7_UT_D5 = condition_IDP13D5Untreated - condition_IDP7D5Untreated,
  P13vsP5_UT_D5 = condition_IDP13D5Untreated - condition_IDP5D5Untreated,
  P7vsP5_T_D5 = condition_IDP7D5Treated - condition_IDP5D5Treated,
  P13vsP7_T_D5 = condition_IDP13D5Treated - condition_IDP7D5Treated,
  P13vsP5_T_D5 = condition_IDP13D5Treated - condition_IDP5D5Treated,

  ## Coefs 25 - 30: Passage x Treatment at each day
  P7vsP5_TvsUT_D3 = (condition_IDP7D3Treated - condition_IDP7D3Untreated) -
    (condition_IDP5D3Treated - 0),
  P13vsP7_TvsUT_D3 = (condition_IDP13D3Treated - condition_IDP13D3Untreated) -
    (condition_IDP7D3Treated - condition_IDP7D3Untreated),
  P13vsP5_TvsUT_D3 = (condition_IDP13D3Treated - condition_IDP13D3Untreated) -
    (condition_IDP5D3Treated - 0),
  P7vsP5_TvsUT_D5 = (condition_IDP7D5Treated - condition_IDP7D5Untreated) -
    (condition_IDP5D5Treated - condition_IDP5D5Untreated),
  P13vsP7_TvsUT_D5 = (condition_IDP13D5Treated - condition_IDP13D5Untreated) -
    (condition_IDP7D5Treated - condition_IDP7D5Untreated),
  P13vsP5_TvsUT_D5 = (condition_IDP13D5Treated - condition_IDP13D5Untreated) -
    (condition_IDP5D5Treated - condition_IDP5D5Untreated),

  ## Coefs 31 - 36: Passage x Day at treatment
  P7vsP5_D5vsD3_UT = (condition_IDP7D5Untreated - condition_IDP5D5Untreated) -
    (condition_IDP7D3Untreated - 0),
  P13vsP7_D5vsD3_UT = (condition_IDP13D5Untreated - condition_IDP7D5Untreated) -
    (condition_IDP13D3Untreated - condition_IDP7D3Untreated),
  P13vsP5_D5vsD3_UT = (condition_IDP13D5Untreated - condition_IDP5D5Untreated) -
    (condition_IDP13D3Untreated - 0),
  P7vsP5_D5vsD3_T = (condition_IDP7D5Treated - condition_IDP5D5Treated) -
    (condition_IDP7D3Treated - condition_IDP5D3Treated),
  P13vsP7_D5vsD3_T = (condition_IDP13D5Treated - condition_IDP7D5Treated) -
    (condition_IDP13D3Treated - condition_IDP7D3Treated),
  P13vsP5_D5vsD3_T = (condition_IDP13D5Treated - condition_IDP5D5Treated) -
    (condition_IDP13D3Treated - condition_IDP5D3Treated),

  ## Coefs 37 - 39: Day x Treatment at Passage
  D5vsD3_TvsUT_P5 = (condition_IDP5D5Treated - condition_IDP5D5Untreated) -
    (condition_IDP5D3Treated - 0),
  D5vsD3_TvsUT_P7 = (condition_IDP7D5Treated - condition_IDP7D5Untreated) -
    (condition_IDP7D3Treated - condition_IDP7D3Untreated),
  D5vsD3_TvsUT_P13 = (condition_IDP13D5Treated - condition_IDP13D5Untreated) -
    (condition_IDP13D3Treated - condition_IDP13D3Untreated),
  levels = design
)

# Define function for converting gene set GMT files to limma::camera format
get_camera_gs <- function(con) {
  camera_gs <- getGmt(con = con) %>%
    geneIds() %>%
    ids2indices(identifiers = rlog_camera$genes$entrezid)

  return(camera_gs)
}

# Obtain gene sets
## List of gene sets to be used for camera
## MSigDB h
## MSigDB C2/CGP
## MSigDB C2/CP
## GO (CC, BP, MF)
list_gmt_camera <- list_gmt %>%
  map(\(x) get_camera_gs(con = x))

# Define camera helper function
# Function runs camera on MArrayLM fit object provided by voomLmFit,
# allowing selection of contrast coefficient and threshold
# Export files to relevant folder
run_camera <- function(
  elist = NULL,
  design = NULL,
  contrasts = NULL,
  genesets = NULL,
  coef = NULL,
  inter.gene.cor = 0.01, # nolint: object_name_linter.
  sort = TRUE
) {
  # Check for errors in input
  ## Check errors in inter-gene correlation
  if (is.null(inter.gene.cor) == TRUE) {
    message(str_c("inter-gene correlation estimated by camera"))
  } else {
    if (
      (is.numeric(inter.gene.cor) == FALSE ||
        abs(inter.gene.cor) > 1)
    ) {
      stop("inter-gene correlation must be between -1 and 1")
    }
    message(str_c("inter-gene correlation set at ", inter.gene.cor))
  }

  ## Check errors in sort option
  if (class(sort) != "logical") {
    stop("sort selection must be TRUE or FALSE")
  }

  # Set up output folder
  path_output <- here::here(
    "output",
    "data_enrichment",
    "camera",
    as.character(colnames(contrasts)[coef])
  )

  message(str_c(
    "Output camera results to",
    path_output,
    sep = " "
  ))

  if (!dir.exists(path_output)) {
    dir.create(path_output, recursive = TRUE)
  }

  # Run camera on list of provided gene sets
  camera_results <- imap(
    list_gmt_camera,
    \(geneset, genesetname) {
      message(str_c(
        "Running",
        genesetname,
        "enrichment for",
        path_output,
        sep = " "
      ))

      camera(
        elist$E,
        geneset,
        design,
        contrasts[, coef],
        sort = sort,
        inter.gene.cor = inter.gene.cor
      )
    }
  )

  # Save data
  saveRDS(
    camera_results,
    file = here::here(
      path_output,
      "camera_results.RDS"
    )
  )

  write.xlsx(
    camera_results,
    file = here::here(
      path_output,
      "camera_results.xlsx"
    ),
    asTable = TRUE,
    rowNames = TRUE
  )

  # Return
  return(camera_results)
}
