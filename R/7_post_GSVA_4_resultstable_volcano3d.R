here::i_am("R/7_post_GSVA_4_resultstable_volcano3d.R")

########################
# Extract results tables for GSVA volcano3D plots
########################

# Import packages
library(conflicted)
library(DESeq2)
library(GSVA)
library(here)
library(magrittr)
library(tidyverse)
library(WGCNA)

# Load data
fit_gsva <- readRDS(
  here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "GSVA_results.RDS"
  )
)

quant_gsva <- readRDS(
  here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "quant_GSVA.RDS"
  )
)

polar_manual_gsva <- readRDS(
  file = here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "polar_manual_gsva.RDS"
  )
)

# Define table specs
## Number of top gene sets
ntop <- 5

# Define coefficients for extracting genes
coef_volcano3d <- list(
  up = list(
    p5 = "5+",
    p7 = "7+",
    p13 = "1+"
  ),
  down = list(
    p5 = "7+1+",
    p7 = "5+1+",
    p13 = "5+7+"
  )
)

coef_toptable <- list(
  p5 = c(p5_p7 = "P7vsP5_UT_D3", p5_p13 = "P13vsP5_UT_D3"),
  p7 = c(p5_p7 = "P7vsP5_UT_D3", p7_p13 = "P13vsP7_UT_D3"),
  p13 = c(p7_p13 = "P13vsP7_UT_D3", p5_p13 = "P13vsP5_UT_D3")
)

# Define collections to consider
collections_interest <- c(
  "h",
  "c2_cgp",
  "c2_cp",
  "GOBP",
  "GOCC",
  "GOMF"
)

# Extract gene set IDs relevant
genesetid_gsva <- purrr::map(polar_manual_gsva[1:6], \(polar_manual) {
  purrr::map(coef_volcano3d, \(x) {
    purrr::map(x, \(coefs) {
      volcano3D::significance_subset(
        polar_manual,
        significance = coefs,
        output = "pvals"
      ) %>%
        rownames()
    })
  })
})

genesetid_wgcna <- purrr::map(coef_volcano3d, \(x) {
  purrr::map(x, \(coefs) {
    polar_manual_gsva[["WGCNA"]]@df$scaled %>%
      dplyr::filter(lab == coefs) %>%
      rownames()
  })
})

genesetid_gsva <- c(genesetid_gsva, list("WGCNA" = genesetid_wgcna))

# Merge fit_gsva lists to get list of
# genesets up/downregulated at each passage compared to the rest
toptable_paired <- purrr::imap(fit_gsva, \(fit, name) {
  purrr::map(coef_toptable, \(contrast) {
    extract_joined_results_limma(
      fit,
      contrast_1 = contrast[[1]],
      contrast_2 = contrast[[2]]
    ) %>%
      dplyr::mutate(collection = name)
  })
})

# Select for gene sets highlighted in 3D volcano plots
toptable_ordered <- purrr::map2(
  toptable_paired,
  genesetid_gsva,
  \(toptable_collection, genesetid_collection) {
    purrr::map(genesetid_collection, \(genesetid) {
      purrr::map2(
        toptable_collection,
        genesetid,
        \(top, genesetid) {
          dplyr::filter(top, geneset %in% genesetid) %>%
            dplyr::arrange(., desc(.[[2]]))
        }
      )
    })
  }
)

# Merge into one table for export to thesis
toptable_out <- toptable_ordered %$%
  list(
    p5 = rbind(
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["down"]][["p5"]],
            P7vsP5_UT_D3 > 0,
            P13vsP5_UT_D3 > 0
          ) %>%
            tidyr::drop_na() %>%
            dplyr::slice_head(n = ntop)
        }
      )),
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["up"]][["p5"]],
            P7vsP5_UT_D3 < 0,
            P13vsP5_UT_D3 < 0
          ) %>%
            tidyr::drop_na() %>%
            dplyr::slice_tail(n = ntop)
        }
      ))
    ),
    p7 = rbind(
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["up"]][["p7"]],
            P7vsP5_UT_D3 > 0,
            P13vsP7_UT_D3 < 0
          ) %>%
            tidyr::drop_na() %>%
            dplyr::slice_head(n = ntop)
        }
      )),
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["down"]][["p7"]],
            P7vsP5_UT_D3 < 0,
            P13vsP7_UT_D3 > 0
          ) %>%
            tidyr::drop_na() %>%
            dplyr::slice_tail(n = ntop)
        }
      ))
    ),
    p13 = rbind(
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["up"]][["p13"]],
            P13vsP7_UT_D3 > 0,
            P13vsP5_UT_D3 > 0
          ) %>%
            tidyr::drop_na() %>%
            dplyr::slice_head(n = ntop)
        }
      )),
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["down"]][["p13"]],
            P13vsP7_UT_D3 < 0,
            P13vsP5_UT_D3 < 0
          ) %>%
            tidyr::drop_na() %>%
            dplyr::slice_tail(n = ntop)
        }
      ))
    )
  )

## Same table but no limits to length
toptable_out_big <- toptable_ordered %$%
  list(
    p5 = rbind(
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["down"]][["p5"]],
            P7vsP5_UT_D3 > 0,
            P13vsP5_UT_D3 > 0
          ) %>%
            tidyr::drop_na()
        }
      )),
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["up"]][["p5"]],
            P7vsP5_UT_D3 < 0,
            P13vsP5_UT_D3 < 0
          ) %>%
            tidyr::drop_na()
        }
      ))
    ),
    p7 = rbind(
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["up"]][["p7"]],
            P7vsP5_UT_D3 > 0,
            P13vsP7_UT_D3 < 0
          ) %>%
            tidyr::drop_na()
        }
      )),
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["down"]][["p7"]],
            P7vsP5_UT_D3 < 0,
            P13vsP7_UT_D3 > 0
          ) %>%
            tidyr::drop_na()
        }
      ))
    ),
    p13 = rbind(
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["up"]][["p13"]],
            P13vsP7_UT_D3 > 0,
            P13vsP5_UT_D3 > 0
          ) %>%
            tidyr::drop_na()
        }
      )),
      dplyr::bind_rows(purrr::map(
        collections_interest,
        \(collection) {
          dplyr::filter(
            .[[collection]][["down"]][["p13"]],
            P13vsP7_UT_D3 < 0,
            P13vsP5_UT_D3 < 0
          ) %>%
            tidyr::drop_na()
        }
      ))
    )
  )


# Export data
write.xlsx(
  x = toptable_out,
  asTable = TRUE,
  file <- here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "toptable_volcano3d.xlsx"
  )
)

write.xlsx(
  x = toptable_out_big,
  asTable = TRUE,
  file <- here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "toptable_volcano3d_big.xlsx"
  )
)
