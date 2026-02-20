here::i_am("R/6_post_camera_04_plot_thesis_grid.R")

########################
# Build grid of barcode plots
########################

# Import packages
library(cowplot)
library(here)
library(tidyverse)

## Between passages - highlight on P+13 pathways
## x-axis: 5v7; 7v13; 5v13
## y-axis: E2F targets, G2M checkpoint, IFNa, IFNg, MYCV1, MYCV2
conditions <- c(
  "P7vsP5_UT_D3",
  "P13vsP7_UT_D3",
  "P13vsP5_UT_D3"
)

pathways_p13 <- c(
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT"
)

## Between passages - highlight on P+7 pathways
## x-axisL 5v7; 7v13; 5v13
## y-axis: NFKB, TGFB, Hypoxia
pathways_p7 <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_TGF_BETA_SIGNALING"
)

# Draw plots
## For P13
### Get file paths
files_p13 <- map(
  conditions,
  \(conditions) {
    map(
      pathways_p13,
      \(pathway) {
        here::here(
          "output",
          "data_enrichment",
          "camera",
          conditions,
          "barcodeplots",
          str_c(pathway, ".png", sep = "")
        )
      }
    )
  }
) %>%
  unlist()

### Draw plots
plots_p13 <- plot_grid(
  plotlist = map(
    files_p13,
    \(x) ggdraw() + draw_image(x)
  ),
  ncol = 2
)

### Export plots
ggsave(
  file = here::here(
    "output",
    "data_enrichment",
    "camera",
    "barcodeplots_P13.png"
  ),
  width = 8,
  height = 9
)

## For P7
### Get file paths
files_p7 <- map(
  conditions,
  \(conditions) {
    map(
      pathways_p7,
      \(pathway) {
        here::here(
          "output",
          "data_enrichment",
          "camera",
          conditions,
          "barcodeplots",
          str_c(pathway, ".png", sep = "")
        )
      }
    )
  }
) %>%
  unlist()

### Draw plots
plots_p7 <- plot_grid(
  plotlist = map(
    files_p7,
    \(x) ggdraw() + draw_image(x)
  ),
  ncol = 2
)

### Export plots
ggsave(
  file = here::here(
    "output",
    "data_enrichment",
    "camera",
    "barcodeplots_P7.png"
  ),
  width = 8,
  height = 9
)
