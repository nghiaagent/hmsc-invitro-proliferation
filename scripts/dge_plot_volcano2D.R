# Load data

source("./scripts/dge_cellpops_as_fixed.R")

# Extract logFC and p-vals for comparisons between passages at D3 and D5

top_P7vsP5_D3 <- topTable(fit_contrasts, coef = 4,
                          number = Inf)

top_P7vsP5_D5 <- topTable(fit_contrasts, coef = 13,
                          number = Inf)

top_P13vsP5_D3 <- topTable(fit_contrasts, coef = 5,
                           number = Inf)

top_P13vsP5_D5 <- topTable(fit_contrasts, coef = 14,
                           number = Inf)

top_P13vsP7_D3 <- topTable(fit_contrasts, coef = 6,
                           number = Inf)

top_P13vsP7_D5 <- topTable(fit_contrasts, coef = 15,
                           number = Inf)


# Extract logFC and p-vals for Trt vs Ctrl

top_treat_P5_D3 <- topTable(fit_contrasts, coef = 1,
                            number = Inf)

top_treat_P5_D5 <- topTable(fit_contrasts, coef = 10,
                            number = Inf)

top_treat_P7_D3 <- topTable(fit_contrasts, coef = 2,
                            number = Inf)

top_treat_P7_D5 <- topTable(fit_contrasts, coef = 11,
                            number = Inf)

top_treat_P13_D3 <- topTable(fit_contrasts, coef = 3,
                             number = Inf)

top_treat_P13_D5 <- topTable(fit_contrasts, coef = 12,
                             number = Inf)

# Draw 2D volcano plot - between growth phases

volcano_P7vsP5_D3 <- EnhancedVolcano::EnhancedVolcano(
  top_P7vsP5_D3,
  lab = top_P7vsP5_D3$GENENAME,
  x = 'logFC',
  y = 'adj.P.Val',
  title = "Phase B - Phase A, after 3 days in culture",
  subtitle = NULL,
  xlim = c(-16.5, 16.5),
  ylim = c(0, 15),
  pCutoff = 0.05,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors = 'black',
  caption = NULL,
  legendDropLevels = F
)


volcano_P7vsP5_D5 <- EnhancedVolcano::EnhancedVolcano(
  top_P7vsP5_D5,
  lab = top_P7vsP5_D5$GENENAME,
  x = 'logFC',
  y = 'adj.P.Val',
  title = "Phase B - Phase A, after 5 days in culture",
  subtitle = NULL,
  xlim = c(-8, 25),
  ylim = c(0, 15),
  pCutoff = 0.05,
  legendPosition = 'right',
  boxedLabels = TRUE,
  drawConnectors = TRUE,
  widthConnectors = 1.0,
  colConnectors = 'black',
  caption = NULL,
  legendDropLevels = F
)


volcano_P13vsP5_D3 <-
  EnhancedVolcano::EnhancedVolcano(
    top_P13vsP5_D3,
    lab = top_P13vsP5_D3$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Phase C - Phase A, after 3 days in culture",
    subtitle = NULL,
    xlim = c(-16.5, 16.5),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )


volcano_P13vsP5_D5 <-
  EnhancedVolcano::EnhancedVolcano(
    top_P13vsP5_D5,
    lab = top_P13vsP5_D5$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Phase C - Phase A, after 5 days in culture",
    subtitle = NULL,
    xlim = c(-8, 25),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

volcano_P13vsP7_D3 <-
  EnhancedVolcano::EnhancedVolcano(
    top_P13vsP7_D3,
    lab = top_P13vsP7_D3$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Phase C - Phase B, after 3 days in culture",
    subtitle = NULL,
    xlim = c(-16.5, 16.5),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )


volcano_P13vsP7_D5 <-
  EnhancedVolcano::EnhancedVolcano(
    top_P13vsP7_D5,
    lab = top_P13vsP7_D5$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Phase C - Phase B, after 5 days in culture",
    subtitle = NULL,
    xlim = c(-8, 25),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )
# Draw 2D volcano plot - between treatments

volcano_treat_P5_D3 <-
  EnhancedVolcano::EnhancedVolcano(
    top_treat_P5_D3,
    lab = top_treat_P5_D3$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Effect of Heparin treatment at Phase A, Day 3",
    subtitle = NULL,
    xlim = c(-9.5, 23.5),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

volcano_treat_P7_D3 <-
  EnhancedVolcano::EnhancedVolcano(
    top_treat_P7_D3,
    lab = top_treat_P7_D3$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Effect of Heparin treatment at Phase B, Day 3",
    subtitle = NULL,
    xlim = c(-9.5, 23.5),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

volcano_treat_P13_D3 <-
  EnhancedVolcano::EnhancedVolcano(
    top_treat_P13_D3,
    lab = top_treat_P13_D3$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Effect of Heparin treatment at Phase C, Day 3",
    subtitle = NULL,
    xlim = c(-9.5, 23.5),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )


volcano_treat_P5_D5 <-
  EnhancedVolcano::EnhancedVolcano(
    top_treat_P5_D5,
    lab = top_treat_P5_D5$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Effect of Heparin treatment at Phase A, Day 5",
    subtitle = NULL,
    xlim = c(-16.5, 16.5),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

volcano_treat_P7_D5 <-
  EnhancedVolcano::EnhancedVolcano(
    top_treat_P7_D5,
    lab = top_treat_P7_D5$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Effect of Heparin treatment at Phase B, Day 5",
    subtitle = NULL,
    xlim = c(-16.5, 16.5),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

volcano_treat_P13_D5 <-
  EnhancedVolcano::EnhancedVolcano(
    top_treat_P13_D5,
    lab = top_treat_P13_D5$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Effect of Heparin treatment at Phase C, Day 5",
    subtitle = NULL,
    xlim = c(-16.5, 16.5),
    ylim = c(0, 15),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )



# Arrange plots in a grid

grid_day3 <-
  plot_grid(
    volcano_P7vsP5_D3  + theme(legend.position = "none"),
    volcano_P13vsP5_D3  + theme(legend.position = "none"),
    volcano_P13vsP7_D3 + theme(legend.position = "bottom"),
    labels = "AUTO",
    label_size = 28,
    ncol = 1,
    nrow = 3
  )

grid_day5 <-
  plot_grid(
    volcano_P7vsP5_D5  + theme(legend.position = "none"),
    volcano_P13vsP5_D5  + theme(legend.position = "none"),
    volcano_P13vsP7_D5 + theme(legend.position = "bottom"),
    labels = "AUTO",
    label_size = 28,
    ncol = 1,
    nrow = 3
  )

grid_treat_day3 <-
  plot_grid(
    volcano_treat_P5_D3  + theme(legend.position = "none"),
    volcano_treat_P7_D3  + theme(legend.position = "none"),
    volcano_treat_P13_D3 + theme(legend.position = "bottom"),
    labels = "AUTO",
    label_size = 28,
    ncol = 1,
    nrow = 3
  )

grid_treat_day5 <-
  plot_grid(
    volcano_treat_P5_D5  + theme(legend.position = "none"),
    volcano_treat_P7_D5  + theme(legend.position = "none"),
    volcano_treat_P13_D5 + theme(legend.position = "bottom"),
    labels = "AUTO",
    label_size = 28,
    ncol = 1,
    nrow = 3
  )

ggsave(
  filename = "./output/plots_volcano/day3.png",
  grid_day3,
  scale = 1.5,
  width = 6,
  height = 10,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "./output/plots_volcano/day5.png",
  grid_day5,
  scale = 1.5,
  width = 6,
  height = 10,
  units = "in",
  dpi = 300
)


ggsave(
  filename = "./output/plots_volcano/treat_day3.png",
  grid_treat_day3,
  scale = 1.5,
  width = 6,
  height = 10,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "./output/plots_volcano/treat_day5.png",
  grid_treat_day5,
  scale = 1.5,
  width = 6,
  height = 10,
  units = "in",
  dpi = 300
)