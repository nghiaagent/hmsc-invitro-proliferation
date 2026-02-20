here::i_am("R/7_post_GSVA_4_resultstable.R")

########################
# Load data
########################

# Import packages
library(DESeq2)
library(here)
library(limma)
library(magrittr)
library(tidyverse)

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

# Set conditions of interest
conditions_interest <- c(
  ## Coefs 1 - 6: Treatment at each timepoint
  "Trt_P5_D3" = 1,
  "Trt_P5_D5" = 2,
  "Trt_P7_D3" = 3,
  "Trt_P7_D5" = 4,
  "Trt_P13_D3" = 5,
  "Trt_P13_D5" = 6,

  ## Coefs 7 - 12: Day at each timepoint x treatment
  "D5vsD3_UT_P5" = 7,
  "D5vsD3_UT_P7" = 8,
  "D5vsD3_UT_P13" = 9,
  "D5vsD3_T_P5" = 10,
  "D5vsD3_T_P7" = 11,
  "D5vsD3_T_P13" = 12,

  ## Coefs 13 - 24: Passage at each day x treatment
  "P7vsP5_UT_D3" = 13,
  "P13vsP7_UT_D3" = 14,
  "P13vsP5_UT_D3" = 15,
  "P7vsP5_T_D3" = 16,
  "P13vsP7_T_D3" = 17,
  "P13vsP5_T_D3" = 18,
  "P7vsP5_UT_D5" = 19,
  "P13vsP7_UT_D5" = 20,
  "P13vsP5_UT_D5" = 21,
  "P7vsP5_T_D5" = 22,
  "P13vsP7_T_D5" = 23,
  "P13vsP5_T_D5" = 24
)

# Create summary tables
results_summary <- map(
  fit_gsva,
  \(x) summary(decideTests(x))
)

# Create nested list containing TopTables
toptables_gsva <- imap(
  fit_gsva,
  \(fit, collection) {
    map(
      conditions_interest,
      \(coef) {
        topTable(
          fit,
          coef = coef,
          number = Inf,
          sort.by = "p"
        ) %>%
          rownames_to_column(var = "set") %>%
          mutate(source = collection)
      }
    )
  },
  .progress = TRUE
)

# Merge nested into one list
toptables_gsva_format <- map(
  conditions_interest,
  \(condition) {
    toptables_gsva %$%
      rbind(
        .[["h"]][[condition]][1:30, ],
        .[["c2_cgp"]][[condition]][1:30, ],
        .[["c2_cp"]][[condition]][1:30, ],
        .[["GOBP"]][[condition]][1:30, ],
        .[["GOCC"]][[condition]][1:30, ],
        .[["GOMF"]][[condition]][1:30, ]
      ) %>%
      dplyr::select(
        set,
        logFC,
        adj.P.Val,
        source
      ) %>%
      dplyr::rename(
        "Gene set name" = "set",
        "LogFC" = "logFC",
        "adj. P-val" = "adj.P.Val",
        "Collection" = "source"
      )
  }
)

# Save data
## Export summary
write.xlsx(
  results_summary,
  file = here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "GSVA_limma_summary.xlsx"
  )
)

## Export top tables
write.xlsx(
  toptables_gsva_format,
  file = here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "GSVA_limma_toptable.xlsx"
  ),
  asTable = TRUE
)
