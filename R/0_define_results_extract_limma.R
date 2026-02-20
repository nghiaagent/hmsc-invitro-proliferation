## Create function to join 2 toptables from a fit object and contrasts
extract_joined_results_limma <- function(
  fit,
  contrast_1,
  contrast_2
) {
  # Format results
  ## Select columns to keep
  vars_to_keep <- c(contrast_1, contrast_2, "adj.P.Val")

  # Construct table
  toptable_anova <- topTable(
    fit,
    coef = NULL,
    number = Inf,
    sort.by = "none"
  ) %>%
    dplyr::select(all_of(vars_to_keep)) %>%
    rownames_to_column(var = "geneset")

  # Return data
  return(toptable_anova)
}
