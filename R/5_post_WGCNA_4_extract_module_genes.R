here::i_am("R/5_post_WGCNA_4_extract_module_genes.R")

########################
# Extract module genes to be analysed further
########################

# Import packages
library(conflicted)
library(data.table)
library(here)
library(tidyverse)
library(WGCNA)

# Load data
## GCN with all samples
gcn <- readRDS(
  file = here::here(
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "GCN_output.RDS"
  )
)

## Add module assignment to gene list
gcn$genes$module <- gcn$net$colors %>%
  WGCNA::labels2colors()

gcn$genes <- gcn$genes %>%
  tibble::rownames_to_column(var = "entrezid_unique")

## Extract gene lists
## entrez for GSVA
## Symbol for String/PPI
gcn_genelists_entrez <- split(
  gcn$genes,
  gcn$genes$module,
  drop = TRUE
) %>%
  purrr::map(
    \(x) {
      dplyr::select(x, entrezid_unique) %>%
        purrr::as_vector() %>%
        unname()
    }
  )

gcn_genelists_symbol <- split(
  gcn$genes,
  gcn$genes$module,
  drop = TRUE
) %>%
  purrr::map(
    \(x) {
      dplyr::select(x, entrezid_unique) %>%
        purrr::as_vector() %>%
        unname()
    }
  )

names(gcn_genelists_entrez) <- stringr::str_c(
  "WGCNA_allsamples_",
  names(gcn_genelists_entrez)
)

names(gcn_genelists_symbol) <- stringr::str_c(
  "WGCNA_allsamples_",
  names(gcn_genelists_symbol)
)

# Define function to write list to GMT file
### Function by Levi Waldron.
### R list object that will be converted to GMT file.  Each element
### should contain a vector of gene names, and the names of the
### elements will used for the gene set names
### Output file name for .gmt file
write_gmt <- function(object, fname) {
  if (class(object) != "list") {
    stop("object should be of class 'list'")
  }

  if (file.exists(fname)) {
    unlink(fname)
  }

  for (i in seq_along(object)) {
    write.table(
      t(c(
        make.names(rep(names(object)[i], 2)),
        object[[i]]
      )),
      sep = "\t",
      quote = FALSE,
      file = fname,
      append = TRUE,
      col.names = FALSE,
      row.names = FALSE
    )
  }
}

# Write gene list in ENTREZID to GMT file
write_gmt(
  gcn_genelists_entrez,
  here::here(
    ".",
    "input",
    "genesets",
    "gcn_sets_WGCNA.gmt"
  )
)

# Write gene list with GENENAME to text files for STRING
purrr::imap(
  list(
    "turquoise_genes" = "WGCNA_allsamples_turquoise",
    "pink_genes" = "WGCNA_allsamples_pink",
    "blue_genes" = "WGCNA_allsamples_blue",
    "red_genes" = "WGCNA_allsamples_red",
    "cyan_genes" = "WGCNA_allsamples_cyan"
  ),
  \(table, name) {
    data.table::fwrite(
      list(gcn_genelists_symbol[[table]]),
      file = here::here(
        "output",
        "data_WGCNA",
        "WGCNA_allsamples",
        str_c(name, ".txt")
      ),
      na = ""
    )
  }
)
