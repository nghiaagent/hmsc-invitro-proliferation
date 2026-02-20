conflicted::conflict_prefer("setdiff", "base")

source("renv/activate.R")

library(tidyverse)
library(BiocParallel)

# Source local function files
## Results formatting
## Results extraction
## Results clipping for plotting with EnhancedVolcano
c(
  "0_define_results_format.R",
  "0_define_results_extract.R",
  "0_define_results_extract_limma.R",
  "0_define_results_clip.R",
  "0_define_GOIs.R",
  "0_define_colours.R",
  "0_define_contrasts.R"
) %>%
  purrr::walk(\(x) {
    message(paste0("Sourcing ", here::here("scripts", x)))
    source(here::here("R", x), echo = FALSE, verbose = FALSE)
  })

# Limit number of cores used due to memory issues on laptops.
BiocParallel::register(
  SnowParam(workers = 1),
  default = TRUE
)

BiocParallel::bpparam()
