# Ideally, this would run after dge_selected_genes for maximum power.
# I choose to run it this way so we can plot genes that are unexpressed.
# Source env_prep.R first.

source("./scripts/dge_plot_volcano3D.R")

# Plot GPCs

plot_GPC1 <- boxplot_trio(polar = polar_manual,
                          value = "GPC1") +
  scale_y_continuous(limits = c(-6,8),
                     n.breaks = 15) +
  theme_bw()


plot_GPC2 <- boxplot_trio(polar = polar_manual,
                          value = "GPC2") +
  scale_y_continuous(limits = c(-6,8),
                     n.breaks = 15) +
  theme_bw()


plot_GPC3 <- boxplot_trio(polar = polar_manual,
                          value = "GPC3") +
  scale_y_continuous(limits = c(-6,8),
                     n.breaks = 15) +
  theme_bw()


plot_GPC4 <- boxplot_trio(polar = polar_manual,
                          value = "GPC4") +
  scale_y_continuous(limits = c(-6,8),
                     n.breaks = 15) +
  theme_bw()


plot_GPC5 <- boxplot_trio(polar = polar_manual,
                          value = "GPC5") +
  scale_y_continuous(limits = c(-6,8),
                     n.breaks = 15) +
  theme_bw()


plot_GPC6 <- boxplot_trio(polar = polar_manual,
                          value = "GPC6") +
  scale_y_continuous(limits = c(-6,8),
                     n.breaks = 15) +
  theme_bw()

plot_grid_GPC <- plot_grid(
  plot_GPC1  + theme(legend.position = "none"),
  plot_GPC2  + theme(legend.position = "none"),
  plot_GPC3  + theme(legend.position = "none"),
  plot_GPC4  + theme(legend.position = "none"),
  plot_GPC5  + theme(legend.position = "none"),
  plot_GPC6  + theme(legend.position = "none"),
  nrow = 2
)

# Plot SDCs

plot_SDC1 <- boxplot_trio(polar = polar_manual,
                          value = "SDC1") +
  scale_y_continuous(limits = c(3.5,10),
                     n.breaks = 6) +
  theme_bw()

plot_SDC2 <- boxplot_trio(polar = polar_manual,
                          value = "SDC2") +
  scale_y_continuous(limits = c(3.5,10),
                     n.breaks = 6) +
  theme_bw()

plot_SDC3 <- boxplot_trio(polar = polar_manual,
                          value = "SDC3") +
  scale_y_continuous(limits = c(3.5,10),
                     n.breaks = 6) +
  theme_bw()

plot_SDC4 <- boxplot_trio(polar = polar_manual,
                          value = "SDC4") +
  scale_y_continuous(limits = c(3.5,10),
                     n.breaks = 6) +
  theme_bw()

plot_grid_SDC <- plot_grid(
  plot_SDC1  + theme(legend.position = "none"),
  plot_SDC2  + theme(legend.position = "none"),
  plot_SDC3  + theme(legend.position = "none"),
  plot_SDC4  + theme(legend.position = "none"),
  nrow = 2
)

# Plot E2Fs

# Export plots

ggsave(filename = "./output/plots_boxplot/GPC.png",
       plot = plot_grid_GPC,
       scale = 1.2,
       width = 10,
       height = 6
       )

ggsave(filename = "./output/plots_boxplot/SDC.png",
       plot = plot_grid_SDC,
       scale = 1.2,
       width = 10,
       height = 6
)
