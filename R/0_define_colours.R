here::i_am("R/0_define_colours.R")

########################
# Define palettes
########################

# Import packages
library(colorspace)
library(conflicted)
library(grDevices)
library(here)
library(tidyverse)
library(vctrs)

# Palette 1: Set of colours for cell lines and timepoints
# hMSC P5, hMSC P7, hMSC P13, RCX, RVM, 1321, SH (7 colours) + treatment variant
# For Figs 1C, 2D, 3G, 5
## Get colours for palettes
## Extract colours from R Okabe-Ito palette in desired order)
## Consistent with Volcano3D plot
### P5: Blue
### P7: Red (Vermillion)
### P13: Green
### RCX: Yellow
### RVM: Orange
### 1321: Purple
### SH: Grey
palette <- c(
  grDevices::palette.colors()[[6]],
  grDevices::palette.colors()[[7]],
  grDevices::palette.colors()[[4]],
  grDevices::palette.colors()[[5]],
  grDevices::palette.colors()[[2]],
  grDevices::palette.colors()[[8]],
  grDevices::palette.colors()[[9]]
)

palette_merge <- vctrs::vec_interleave(
  colorspace::lighten(palette, 0.2),
  colorspace::darken(palette, 0.2),
)

# Palette 2: Set of colours for heatmap
# For Figs 1A
## Scale 1: Cell population
## Scale 2: Passage
## Scale 3: Day
## Scale 4: Treatment
palette_heatmap <- list(
  `Cell population` = c(
    "hMSC-20176" = grDevices::palette.colors()[1],
    "hMSC-21558" = grDevices::palette.colors()[9]
  ),
  `Passage` = c(
    "P5" = palette_merge[[1]],
    "P7" = palette_merge[[3]],
    "P13" = palette_merge[[5]]
  ),
  `Day` = c(
    "D3" = grDevices::hcl.colors(5, palette = "Purples 3")[3],
    "D5" = grDevices::hcl.colors(5, palette = "Purples 3")[1]
  ),
  `Treatment` = c(
    "Control" = grDevices::hcl.colors(5, palette = "Reds 3")[3],
    "Heparin" = grDevices::hcl.colors(5, palette = "Reds 3")[1]
  ),
  `Batch` = c(
    "IanAmpliseqTranscriptome550Chip1434" = grDevices::palette.colors(
      palette = "Polychrome 36"
    )[[3]],
    "IanAmpliseqTranscriptome550Chip2435" = grDevices::palette.colors(
      palette = "Polychrome 36"
    )[[4]],
    "IanAmpliseqtranscriptome550Chip3437" = grDevices::palette.colors(
      palette = "Polychrome 36"
    )[[5]],
    "R20171107152458" = grDevices::palette.colors(palette = "Polychrome 36")[[
      6
    ]],
    "R20171202120807" = grDevices::palette.colors(palette = "Polychrome 36")[[
      7
    ]],
    "R20171204095805" = grDevices::palette.colors(palette = "Polychrome 36")[[
      8
    ]],
    "R20171204123655" = grDevices::palette.colors(palette = "Polychrome 36")[[
      9
    ]],
    "R20171205085758" = grDevices::palette.colors(palette = "Polychrome 36")[[
      10
    ]],
    "R20171205113620" = grDevices::palette.colors(palette = "Polychrome 36")[[
      11
    ]],
    "R20180131074914" = grDevices::palette.colors(palette = "Polychrome 36")[[
      12
    ]]
  )
)

# Palette 3: For Volcano3D
palette_volcano3d <- c(
  grDevices::palette.colors()[[9]],
  grDevices::palette.colors()[[7]],
  grDevices::palette.colors()[[5]],
  grDevices::palette.colors()[[4]],
  grDevices::palette.colors()[[3]],
  grDevices::palette.colors()[[6]],
  grDevices::palette.colors()[[8]]
)

# Palette 3: For Fig S1 (PCA)
palette_pca <- list(
  "run_date" = grDevices::palette.colors(palette = "Polychrome 36"),
  "cell_line" = c(
    "hMSC-20176" = grDevices::palette.colors()[1],
    "hMSC-21558" = grDevices::palette.colors()[9]
  ),
  "timepoint_ID" = c(
    P5D3 = palette_merge[[1]],
    P5D5 = palette_merge[[2]],
    P7D3 = palette_merge[[3]],
    P7D5 = palette_merge[[4]],
    P13D3 = palette_merge[[5]],
    P13D5 = palette_merge[[6]]
  ),
  "Treatment" = c(
    "Untreated" = grDevices::hcl.colors(5, palette = "Reds 3")[3],
    "Treated" = grDevices::hcl.colors(5, palette = "Reds 3")[1]
  )
)
# Palette 4: For Fig 4 to make quadrant genes more distinct
palette_quadrant <- c(
  "Control only" = palette()[[4]],
  "Heparin only" = palette()[[2]],
  "Both" = palette()[[8]] %>% colorspace::lighten(0.4),
  "ns" = palette()[[8]] %>% colorspace::darken(0.4)
)

# Palette 5: For WGCNA
palette_wgcna <- c(
  "D3" = grDevices::hcl.colors(5, palette = "Purples 3")[3],
  "D5" = grDevices::hcl.colors(5, palette = "Purples 3")[1],
  "hMSC-20176" = grDevices::palette.colors()[1],
  "hMSC-21558" = grDevices::palette.colors()[9],
  "IanAmpliseqTranscriptome550Chip1434" = grDevices::palette.colors(
    palette = "Polychrome 36"
  )[[3]],
  "IanAmpliseqTranscriptome550Chip2435" = grDevices::palette.colors(
    palette = "Polychrome 36"
  )[[4]],
  "IanAmpliseqtranscriptome550Chip3437" = grDevices::palette.colors(
    palette = "Polychrome 36"
  )[[5]],
  "P13" = palette_merge[[5]],
  "P5" = palette_merge[[1]],
  "P7" = palette_merge[[3]],
  "R20171107152458" = grDevices::palette.colors(palette = "Polychrome 36")[[6]],
  "R20171202120807" = grDevices::palette.colors(palette = "Polychrome 36")[[7]],
  "R20171204095805" = grDevices::palette.colors(palette = "Polychrome 36")[[8]],
  "R20171204123655" = grDevices::palette.colors(palette = "Polychrome 36")[[9]],
  "R20171205085758" = grDevices::palette.colors(palette = "Polychrome 36")[[
    10
  ]],
  "R20171205113620" = grDevices::palette.colors(palette = "Polychrome 36")[[
    11
  ]],
  "R20180131074914" = grDevices::palette.colors(palette = "Polychrome 36")[[
    12
  ]],
  "Treated" = grDevices::hcl.colors(5, palette = "Reds 3")[1],
  "Untreated" = grDevices::hcl.colors(5, palette = "Reds 3")[3]
)
