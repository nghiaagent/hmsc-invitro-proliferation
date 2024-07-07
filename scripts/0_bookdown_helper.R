# Define function for book render

## r Define XY plot methods
# Use glimmaXY to plot only significant gene sets
# Select only significant gene sets
# Generate new signif codes: -1 = x only, 0 = both, 1 = y only


plot_XY <- function(fit,
                    coef_x,
                    coef_y,
                    xlab = "x",
                    ylab = "y") {
  # Add annotations to fit obj
  
  
  # Select only signif gene sets
  
  status <-
    abs(decideTests(fit)[, coef_x]) + abs(decideTests(fit)[, coef_y])
  
  fit_signif <- fit[status > 0, ]
  
  top_x <- topTable(fit_signif,
                    coef = coef_x,
                    number = Inf,
                    sort.by = "none")
  
  top_y <- topTable(fit_signif,
                    coef = coef_y,
                    number = Inf,
                    sort.by = "none")
  
  # Generate signif codes to show overlap
  # How code works
  #    x  DOWN | ns | UP
  #  y
  # DOWN   -4  | -3 | -2
  # ns     -1  |  0 | 1
  # UP     2  |  3 | 4
  # Swap -1 and 1 for  -1 (x only)
  # Swap -3 and 3 for 1 (y only)
  # Swap -4, -2, 2, 4 for 0 (both x and y)
  
  status_signif <-
    (decideTests(fit_signif)[, coef_x] + 3 * decideTests(fit_signif)[, coef_y]) %>%
    case_match(c(-1, 1) ~ -1,
               c(-3, 3) ~ 1,
               c(-4, -2, 2, 4) ~ 0)
  
  
  vec_x <- topTable(fit_signif,
                    coef = coef_x,
                    number = Inf,
                    sort.by = "none")$logFC %>%
    as.data.frame()
  
  rownames(vec_x) <- rownames(fit_signif$coefficients)
  
  vec_y <- topTable(fit_signif,
                    coef = coef_y,
                    number = Inf,
                    sort.by = "none")$logFC %>%
    as.data.frame()
  
  rownames(vec_y) <- rownames(fit_signif$coefficients)
  
  glimmaXY(
    x = vec_x,
    y = vec_y,
    xlab = xlab,
    ylab = ylab,
    status = status_signif,
    anno = fit_signif$genes,
    status.cols = c("#1052bd",
                    "#9951A4",
                    "#cc212f"),
    main = "Blue - WT only, Purple - Both, Red - KO only",
    width = 750,
    height = 900
  )
}

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