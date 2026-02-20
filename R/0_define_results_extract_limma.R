here::i_am("R/0_define_results_extract_limma.R")

########################
# Define function to join 2 toptables from a fit object and contrasts
########################

# Import packages
library(conflicted)
library(DESeq2)
library(here)
library(limma)
library(tidyverse)

# Define function
extract_joined_results_limma <- function(
  fit,
  contrast_1,
  contrast_2
) {
  # Format results
  ## Select columns to keep
  vars_to_keep <- c(contrast_1, contrast_2, "adj.P.Val")

  # Construct table
  toptable_anova <- limma::topTable(
    fit,
    coef = NULL,
    number = Inf,
    sort.by = "none"
  ) %>%
    dplyr::select(dplyr::all_of(vars_to_keep)) %>%
    tibble::rownames_to_column(var = "geneset")

  # Return data
  return(toptable_anova)
}
