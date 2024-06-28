# Extract logFC and p-vals for comparisons between passages at D3 and D5

top_P7vsP5_D3 <- topTable(fit_contrasts, coef = 13,
                          number = Inf)

top_P7vsP5_D5 <- topTable(fit_contrasts, coef = 16,
                          number = Inf)

top_P13vsP5_D3 <- topTable(fit_contrasts, coef = 15,
                           number = Inf)

top_P13vsP5_D5 <- topTable(fit_contrasts, coef = 18,
                           number = Inf)

top_P13vsP7_D3 <- topTable(fit_contrasts, coef = 14,
                           number = Inf)

top_P13vsP7_D5 <- topTable(fit_contrasts, coef = 17,
                           number = Inf)

# Extract logFC and p-vals for Trt vs Ctrl

top_treat_P5_D3 <- topTable(fit_contrasts, coef = 1,
                            number = Inf)

top_treat_P5_D5 <- topTable(fit_contrasts, coef = 2,
                            number = Inf)

top_treat_P7_D3 <- topTable(fit_contrasts, coef = 3,
                            number = Inf)

top_treat_P7_D5 <- topTable(fit_contrasts, coef = 4,
                            number = Inf)

top_treat_P13_D3 <- topTable(fit_contrasts, coef = 5,
                             number = Inf)

top_treat_P13_D5 <- topTable(fit_contrasts, coef = 6,
                             number = Inf)

# Extract logFC and pvals for interaction term

top_interaction_TvsUT_P7vsP5_D3 <-
  topTable(fit_contrasts, coef = 25,
           number = Inf)

top_interaction_TvsUT_P13vsP5_D3 <-
  topTable(fit_contrasts, coef = 27,
           number = Inf)

top_interaction_TvsUT_P13vsP7_D3 <-
  topTable(fit_contrasts, coef = 26,
           number = Inf)

top_interaction_TvsUT_P7vsP5_D5 <-
  topTable(fit_contrasts, coef = 28,
           number = Inf)

top_interaction_TvsUT_P13vsP5_D5 <-
  topTable(fit_contrasts, coef = 30,
           number = Inf)

top_interaction_TvsUT_P13vsP7_D5 <-
  topTable(fit_contrasts, coef = 29,
           number = Inf)

# Draw 2D volcano plot - between growth phases

volcano_P7vsP5_D3 <- EnhancedVolcano::EnhancedVolcano(
  top_P7vsP5_D3,
  lab = top_P7vsP5_D3$GENENAME,
  x = 'logFC',
  y = 'adj.P.Val',
  title = "Phase B - Phase A, after 3 days in culture",
  subtitle = NULL,
  xlim = c(-12, 14),
  ylim = c(0, 17),
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
  xlim = c(-12, 28),
  ylim = c(0, 14),
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
    xlim = c(-12, 14),
    ylim = c(0, 17),
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
    xlim = c(-12, 28),
    ylim = c(0, 14),
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
    xlim = c(-12, 14),
    ylim = c(0, 17),
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
    xlim = c(-12, 28),
    ylim = c(0, 14),
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
    xlim = c(-10, 10),
    ylim = c(0, 3.5),
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
    xlim = c(-10, 10),
    ylim = c(0, 3.5),
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
    xlim = c(-10, 10),
    ylim = c(0, 3.5),
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
    xlim = c(-15, 7),
    ylim = c(0, 5),
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
    xlim = c(-15, 7),
    ylim = c(0, 5),
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
    xlim = c(-15, 7),
    ylim = c(0, 5),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

# Draw 2D volcano plot - interaction term
## Day 3

volcano_interaction_TvsUT_P7vsP5_D3 <-
  EnhancedVolcano::EnhancedVolcano(
    top_interaction_TvsUT_P7vsP5_D3,
    lab = top_interaction_TvsUT_P7vsP5_D3$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Change in effect of Heparin treatment between Phase A & B, Day 3",
    subtitle = NULL,
    xlim = c(-25, 12),
    ylim = c(0, 5),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

volcano_interaction_TvsUT_P13vsP5_D3 <-
  EnhancedVolcano::EnhancedVolcano(
    top_interaction_TvsUT_P13vsP5_D3,
    lab = top_interaction_TvsUT_P13vsP5_D3$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Change in effect of Heparin treatment between Phase A & C, Day 3",
    subtitle = NULL,
    xlim = c(-25, 12),
    ylim = c(0, 5),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

volcano_interaction_TvsUT_P13vsP7_D3 <-
  EnhancedVolcano::EnhancedVolcano(
    top_interaction_TvsUT_P13vsP7_D3,
    lab = top_interaction_TvsUT_P13vsP7_D3$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Change in effect of Heparin treatment between Phase A & C, Day 3",
    subtitle = NULL,
    xlim = c(-25, 12),
    ylim = c(0, 5),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

## Day 5

volcano_interaction_TvsUT_P7vsP5_D5 <-
  EnhancedVolcano::EnhancedVolcano(
    top_interaction_TvsUT_P7vsP5_D5,
    lab = top_interaction_TvsUT_P7vsP5_D5$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Change in effect of Heparin treatment between Phase A & B, Day 5",
    subtitle = NULL,
    xlim = c(-15, 15),
    ylim = c(0, 5),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

volcano_interaction_TvsUT_P13vsP5_D5 <-
  EnhancedVolcano::EnhancedVolcano(
    top_interaction_TvsUT_P13vsP5_D5,
    lab = top_interaction_TvsUT_P13vsP5_D5$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Change in effect of Heparin treatment between Phase A & C, Day 5",
    subtitle = NULL,
    xlim = c(-15, 15),
    ylim = c(0, 5),
    pCutoff = 0.05,
    legendPosition = 'right',
    boxedLabels = TRUE,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    caption = NULL,
    legendDropLevels = F
  )

volcano_interaction_TvsUT_P13vsP7_D5 <-
  EnhancedVolcano::EnhancedVolcano(
    top_interaction_TvsUT_P13vsP7_D5,
    lab = top_interaction_TvsUT_P13vsP7_D5$GENENAME,
    x = 'logFC',
    y = 'adj.P.Val',
    title = "Change in effect of Heparin treatment between Phase A & C, Day 5",
    subtitle = NULL,
    xlim = c(-15, 15),
    ylim = c(0, 5),
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
    volcano_P13vsP7_D3 + theme(legend.position = "none"),
    labels = "AUTO",
    label_size = 28,
    ncol = 1,
    nrow = 3
  )

grid_day5 <-
  plot_grid(
    volcano_P7vsP5_D5  + theme(legend.position = "none"),
    volcano_P13vsP5_D5  + theme(legend.position = "none"),
    volcano_P13vsP7_D5 + theme(legend.position = "none"),
    labels = "AUTO",
    label_size = 28,
    ncol = 1,
    nrow = 3
  )

grid_treat_day3 <-
  plot_grid(
    volcano_treat_P5_D3  + theme(legend.position = "none"),
    volcano_treat_P7_D3  + theme(legend.position = "none"),
    volcano_treat_P13_D3 + theme(legend.position = "none"),
    labels = "AUTO",
    label_size = 28,
    ncol = 1,
    nrow = 3
  )

grid_treat_day5 <-
  plot_grid(
    volcano_treat_P5_D5  + theme(legend.position = "none"),
    volcano_treat_P7_D5  + theme(legend.position = "none"),
    volcano_treat_P13_D5 + theme(legend.position = "none"),
    labels = "AUTO",
    label_size = 28,
    ncol = 1,
    nrow = 3
  )

grid_interaction_day3 <- 
  plot_grid(
    volcano_interaction_TvsUT_P7vsP5_D3  + theme(legend.position = "none"),
    volcano_interaction_TvsUT_P13vsP5_D3  + theme(legend.position = "none"),
    volcano_interaction_TvsUT_P13vsP7_D3 + theme(legend.position = "none"),
    labels = "AUTO",
    label_size = 28,
    ncol = 1,
    nrow = 3
  )

grid_interaction_day5 <- 
  plot_grid(
    volcano_interaction_TvsUT_P7vsP5_D5  + theme(legend.position = "none"),
    volcano_interaction_TvsUT_P13vsP5_D5  + theme(legend.position = "none"),
    volcano_interaction_TvsUT_P13vsP7_D5 + theme(legend.position = "none"),
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

ggsave(
  filename = "./output/plots_volcano/interaction_day3.png",
  grid_interaction_day3,
  scale = 1.5,
  width = 6,
  height = 10,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "./output/plots_volcano/interaction_day5.png",
  grid_interaction_day5,
  scale = 1.5,
  width = 6,
  height = 10,
  units = "in",
  dpi = 300
)

## Make fig for paper

grid_paper_day3 <-
  plot_grid(
    ggdraw() + draw_text("P7vsP5", 
                         angle = 90,
                         fontface = "bold",
                         size = 25),
    volcano_P7vsP5_D3  + 
      theme(legend.position = "right",
            plot.title = element_blank()),
    ggdraw() + draw_text("P13vsP7", 
                         angle = 90,
                         fontface = "bold",
                         size = 25),
    volcano_P13vsP5_D3  + 
      theme(legend.position = "right",
            plot.title = element_blank()),
    ggdraw() + draw_text("P13vsP5", 
                         angle = 90,
                         fontface = "bold",
                         size = 25),
    volcano_P13vsP7_D3 + 
      theme(legend.position = "right",
            plot.title = element_blank()),
    labels = NULL,
    label_size = 28,
    ncol = 2,
    nrow = 3,
    rel_widths = c(0.05,1)
  )

ggsave(
  filename = "./output/plots_volcano/day3_paperver.png",
  grid_paper_day3,
  scale = 0.7,
  width = 18,
  height = 15,
  units = "in",
  dpi = 300
)

## Make fig for paper - between treatments

grid_paper_treat <-
  plot_grid(
    ggdraw() + draw_text("P5D3", 
                         angle = 90,
                         fontface = "bold",
                         size = 25),
    volcano_treat_P5_D3  + 
      theme(legend.position = "right",
            plot.title = element_blank()),
    ggdraw() + draw_text("P7D3", 
                         angle = 90,
                         fontface = "bold",
                         size = 25),
    volcano_treat_P7_D3  + 
      theme(legend.position = "right",
            plot.title = element_blank()),
    ggdraw() + draw_text("P13D3", 
                         angle = 90,
                         fontface = "bold",
                         size = 25),
    volcano_treat_P13_D3 + 
      theme(legend.position = "right",
            plot.title = element_blank()),
    labels = NULL,
    label_size = 28,
    ncol = 2,
    nrow = 3,
    rel_widths = c(0.05,1)
  )

ggsave(
  filename = "./output/plots_volcano/treat_paperver.png",
  grid_paper_treat,
  scale = 0.7,
  width = 18,
  height = 15,
  units = "in",
  dpi = 300
)

### Extra fig for legend

ggsave(
  filename = "./output/plots_volcano/fig_for_legend.png",
  volcano_P7vsP5_D3  + theme(legend.position = "right"),
  scale = 1.5,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)