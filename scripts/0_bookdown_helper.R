# Define function for book render
## Draw camera table

draw_camera_tbl <- function(table) {
  datatable(
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
}