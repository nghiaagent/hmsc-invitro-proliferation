here::i_am("R/4_dge_table_topgenes_volcano3D.R")

########################
# Extract genes from 3D volcano plot
# To put in a table for thesis
########################

# Import packages
library(DESeq2)
library(here)
library(openxlsx)
library(tidyverse)

# Load data
quant_deseq2 <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_batchcor.RDS"
  )
)

results_lrt <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "quant_deseq2_LRT.RDS"
  )
) %>%
  results(alpha = 0.05)

results_lfcshrink <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "post_DGE",
    "results_deseq2_lfcshrink.RDS"
  )
)

polar_manual <- readRDS(
  file = here::here(
    "output",
    "plots_volcano3D",
    "polar_manual.RDS"
  )
)

# Define table specs
ntop <- 20

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

coef_results <- list(
  p5 = c(p5_p7 = 13, p5_p13 = 15),
  p7 = c(p5_p7 = 13, p7_p13 = 14),
  p13 = c(p7_p13 = 14, p5_p13 = 15)
)

# Extract genes based on their direction of regulation
# To find genes uniquely up/downregulated at each passage
# And find genes progressively up/downregulated throughout passages
geneid_volcano3d <- map(
  coef_volcano3d,
  \(x) {
    map(
      x,
      \(coefs) {
        significance_subset(
          polar_manual,
          significance = coefs,
          output = "pvals"
        ) %>%
          rownames()
      }
    )
  }
)

# Get list of genes up/downregulated at each passage compared to the rest
## Single tables
results_volcano3d_single <- map(
  coef_results,
  \(x) {
    extract_joined_results(
      results_1 = results_lfcshrink[[x[[1]]]],
      results_2 = results_lfcshrink[[x[[2]]]],
      results_lrt = results_lrt,
      name_1 = names(x)[[1]],
      name_2 = names(x)[[2]],
      quant_deseq2
    ) %>%
      dplyr::select(c(1, 2, 3, 4, 5, 7)) %>%
      dplyr::rename("adj. P-val" = padj_lrt)
  }
)

results_volcano3d_single <- map(
  geneid_volcano3d,
  \(geneid) {
    map2(
      results_volcano3d_single,
      geneid,
      \(res, id) {
        filter(res, `ENSEMBL ID` %in% id) %>%
          dplyr::arrange(., desc(.[[5]]))
      }
    )
  }
)

## Paired tables
results_volcano3d_pairs <- list(
  p5 = rbind(
    filter(
      results_volcano3d_single[["down"]][["p5"]],
      log2FoldChange_p5_p7 > 0,
      log2FoldChange_p5_p13 > 0
    ),
    filter(
      results_volcano3d_single[["up"]][["p5"]],
      log2FoldChange_p5_p7 < 0,
      log2FoldChange_p5_p13 < 0
    )
  ),
  p7 = rbind(
    filter(
      results_volcano3d_single[["up"]][["p7"]],
      log2FoldChange_p5_p7 > 0,
      log2FoldChange_p7_p13 < 0
    ),
    filter(
      results_volcano3d_single[["down"]][["p7"]],
      log2FoldChange_p5_p7 < 0,
      log2FoldChange_p7_p13 > 0
    )
  ),
  p13 = rbind(
    filter(
      results_volcano3d_single[["up"]][["p13"]],
      log2FoldChange_p7_p13 > 0,
      log2FoldChange_p5_p13 > 0
    ),
    filter(
      results_volcano3d_single[["down"]][["p13"]],
      log2FoldChange_p7_p13 < 0,
      log2FoldChange_p5_p13 < 0
    )
  )
)

results_volcano3d_pairs_top <- list(
  p5 = rbind(
    filter(
      results_volcano3d_single[["down"]][["p5"]],
      log2FoldChange_p5_p7 > 0,
      log2FoldChange_p5_p13 > 0
    ) %>%
      slice_head(n = ntop),
    filter(
      results_volcano3d_single[["up"]][["p5"]],
      log2FoldChange_p5_p7 < 0,
      log2FoldChange_p5_p13 < 0
    ) %>%
      slice_tail(n = ntop)
  ),
  p7 = rbind(
    filter(
      results_volcano3d_single[["up"]][["p7"]],
      log2FoldChange_p5_p7 > 0,
      log2FoldChange_p7_p13 < 0
    ) %>%
      slice_head(n = ntop),
    filter(
      results_volcano3d_single[["down"]][["p7"]],
      log2FoldChange_p5_p7 < 0,
      log2FoldChange_p7_p13 > 0
    ) %>%
      slice_tail(n = ntop)
  ),
  p13 = rbind(
    filter(
      results_volcano3d_single[["up"]][["p13"]],
      log2FoldChange_p7_p13 > 0,
      log2FoldChange_p5_p13 > 0
    ) %>%
      slice_head(n = ntop),
    filter(
      results_volcano3d_single[["down"]][["p13"]],
      log2FoldChange_p7_p13 < 0,
      log2FoldChange_p5_p13 < 0
    ) %>%
      slice_tail(n = ntop)
  )
)

# Export data
write.xlsx(
  x = results_volcano3d_pairs,
  file = here::here(
    "output",
    "data_expression",
    "genes_volcano3d_all.xlsx"
  ),
  asTable = TRUE
)

write.xlsx(
  x = results_volcano3d_pairs_top,
  file = here::here(
    "output",
    "data_expression",
    "genes_volcano3d_top.xlsx"
  ),
  asTable = TRUE
)
