here::i_am("R/0_bookdown_helper.R")

########################
# Define functions for book render
########################

# Import packages
library(DT)
library(here)
library(tidyverse)

# Define how to draw tables for camera results
draw_camera_tbl <- function(table) {
  out <- DT::datatable(
    table,
    options = list(
      columnDefs = list(list(
        targets = 0,
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 30 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
          "}"
        )
      ))
    )
  ) %>%
    # Round p-values to 5 decimal places
    formatRound(c(3, 4), 5)

  return(out)
}
