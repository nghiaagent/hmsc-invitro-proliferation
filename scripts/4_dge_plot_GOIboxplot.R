# Ideally, this would run after dge_selected_genes for maximum power.
# I choose to run it this way so we can plot genes that are unexpressed.
# Source env_prep.R first.

# Select data for conditions

conditions_interest <- c("P5D3Untreated",
                         "P7D3Untreated",
                         "P13D3Untreated")

quant_DGE_small <-
  quant_DGE_voom[, quant_DGE_voom$targets$condition_ID %in% conditions_interest]

plot_GOI_passage <- function(gene_sel) {
  
  gene_sel_name <- mapIds(org.Hs.eg.db,
                                    keys = gene_sel,
                                    column = "SYMBOL",
                                    keytype = "ENTREZID")

quant_sel <-
  quant_DGE_small[quant_DGE_small$genes$ENTREZID == as.integer(gene_sel), ]

df_plot <- quant_sel$targets %>%
  cbind(as.vector(quant_sel$E)) %>%
  rename("logCPM" = "as.vector(quant_sel$E)")

ggplot(df_plot,
       aes(x = Passage,
           y = logCPM,
           fill = Passage)) +
  geom_boxplot() +
  geom_jitter() +
  theme_minimal_hgrid() +
  ggtitle(label = gene_sel_name)

}

plots_inhouse_passage <- map(geneIds_inhouse[geneIds_inhouse %in% as.character(quant_DGE_small$genes$ENTREZID)],
                     plot_GOI_passage)

plots_Wnt_passage <- map(geneIds_Wnt[geneIds_Wnt %in% as.character(quant_DGE_small$genes$ENTREZID)],
                 plot_GOI_passage)

plots_HSPGs_passage <- map(geneIds_HSPGs[geneIds_HSPGs %in% as.character(quant_DGE_small$genes$ENTREZID)],
                   plot_GOI_passage)

plots_rt2array_passage <- map(geneIds_rt2array[geneIds_rt2array %in% as.character(quant_DGE_small$genes$ENTREZID)],
                           plot_GOI_passage)

for (i in 1:length(plots_inhouse_passage)) {
  
  message(str_c("drawing plot for ",
                names(plots_inhouse_passage)[[i]],
                sep = ""))
  
  ggsave(
    filename = str_c(names(plots_inhouse_passage)[[i]],
                           ".png",
                           sep = ""),
    plot = plots_inhouse_passage[[i]],
    path = file.path(
      '.',
      'output',
      'plots_boxplot_logCPM',
      'main'
    ),
    scale = 0.7,
    width = 8,
    height = 8,
    units = "in",
    dpi = 144
  )
}

for (i in 1:length(plots_Wnt_passage)) {
  
  message(str_c("drawing plot for ",
                names(plots_Wnt_passage)[[i]],
                sep = ""))
  
  ggsave(
    filename = str_c(names(plots_Wnt_passage)[[i]],
                     ".png",
                     sep = ""),
    plot = plots_Wnt_passage[[i]],
    path = file.path(
      '.',
      'output',
      'plots_boxplot_logCPM',
      'WNT'
    ),
    scale = 0.7,
    width = 8,
    height = 8,
    units = "in",
    dpi = 144
  )
}

for (i in 1:length(plots_HSPGs_passage)) {
  
  message(str_c("drawing plot for ",
                names(plots_HSPGs_passage)[[i]],
                sep = ""))
  
  ggsave(
    filename = str_c(names(plots_HSPGs_passage)[[i]],
                     ".png",
                     sep = ""),
    plot = plots_HSPGs_passage[[i]],
    path = file.path(
      '.',
      'output',
      'plots_boxplot_logCPM',
      'HSPG_extra'
    ),
    scale = 0.7,
    width = 8,
    height = 8,
    units = "in",
    dpi = 144
  )
}

for (i in 1:length(plots_rt2array_passage)) {
  
  message(str_c("drawing plot for ",
                names(plots_rt2array_passage)[[i]],
                sep = ""))
  
  ggsave(
    filename = str_c(names(plots_rt2array_passage)[[i]],
                     ".png",
                     sep = ""),
    plot = plots_rt2array_passage[[i]],
    path = file.path(
      '.',
      'output',
      'plots_boxplot_logCPM',
      'rt2array'
    ),
    scale = 0.7,
    width = 8,
    height = 8,
    units = "in",
    dpi = 144
  )
}