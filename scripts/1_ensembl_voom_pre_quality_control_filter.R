### Don't source this file by itself; call in from another file after running env_prep.R
### Load dataset 

quant_cDNA_DGE <- readRDS(file = "./output/data_expression/pre_DGE/quant_cDNA_DGE.RDS")

# Library size barplot. Library size is the total number of counts across all genes per sample

## Generate table for this purpose

quant_libsize <- quant_cDNA_DGE$samples %>%
  dplyr::arrange(condition_ID) %>%
  dplyr::select(c("lib.size", "ID", "timepoint_ID", "condition_ID")) %>%
  relocate(ID)

## Barplot

plot_libsize <- ggplot(quant_libsize,
                       aes(
                         x = ID,
                         y = lib.size / 1e6,
                         colour = condition_ID,
                         fill = condition_ID
                       )) +
  geom_bar(stat = "identity",
           alpha = 0.6) +
  labs(x = "",
       y = "Library size (million)",
       colour = "Condition",
       fill = "Condition") +
  scale_y_continuous(breaks = seq(0, 40, 5)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )

ggsave(
  "./output/plots_QC/Library size.png",
  plot = plot_libsize,
  width = 12,
  height = 6,
  scale = 0.7
)

# Get average count for each gene in all samples, sorted from highest to lowest

df_stats_cDNA <- rowMeans(quant_cDNA_DGE$counts) %>%
  as_tibble() %>%
  mutate(ENSEMBL = rownames(quant_cDNA_DGE$counts)) %>%
  relocate(ENSEMBL) %>%
  dplyr::rename("average_counts" = "value") %>%
  mutate(ratio_total_counts = average_counts/sum(average_counts)) %>%
  left_join(quant_cDNA_DGE$genes, by = join_by(ENSEMBL == GENEID)) %>%
  arrange(desc(ratio_total_counts))

## Write list of top 50 genes

head(df_stats_cDNA, 50) %>%
  write_csv("./output/QC/Top_20_Expressed_Genes_cDNA.csv")

## Filter genes with low expression using edgeR function

expr_cutoff <- 40

keep_cDNA_edgeRfiter <- filterByExpr(quant_cDNA_DGE,
                                     group = "condition_ID",
                                     min.count = expr_cutoff,
                                     min.total.count = 1.5*expr_cutoff)

quant_cDNA_DGE_filter <- quant_cDNA_DGE[keep_cDNA_edgeRfiter, , keep.lib.sizes=FALSE] %>%
  calcNormFactors()

# Per sample distribution; before and after adding TMM scaling factor

lcpm <- list("Before normalisation" = quant_cDNA_DGE[keep_cDNA_edgeRfiter, ,
                                                     keep.lib.sizes=FALSE],
             "After normalisation" = quant_cDNA_DGE_filter) %>%
  map2(
    .x = .,
    .y = names(.),
    \ (x, y) {
    cpm(x, log = TRUE) %>%
      as_tibble() %>%
      mutate(GENEID = rownames(quant_cDNA_DGE_filter)) %>%
      pivot_longer(
        cols = !GENEID,
        names_to = "ID",
        values_to = "lcpm"
      ) %>%
      mutate(type = y)
      }
  ) %>%
  bind_rows() %>%
  mutate(type = factor(type, levels = c(
    "Before normalisation",
    "After normalisation"
  ))) %>%
  left_join(y = select(quant_cDNA_DGE_filter$samples, c("ID", "condition_ID")),
            by = join_by(ID == ID)) %>%
  mutate(ID = factor(ID, levels = levels(quant_cDNA_DGE_filter$samples$ID)))

plot_distribution_preTMM <- ggplot(lcpm,
                                   aes(x = ID,
                                       y = lcpm,
                                       colour = condition_ID,
                                       fill = condition_ID)) +
  geom_boxplot(alpha = 0.6) +
  labs(x = "",
       y = "log-counts per million",
       colour = "Condition",
       fill = "Condition") +
  scale_y_continuous(breaks = seq(0, 40, 5)) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_blank(),
    axis.ticks.x=element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  ) +
  facet_wrap( ~ type,
              nrow = 2)

ggsave(
  "./output/plots_QC/distribution.png",
  plot = plot_distribution_preTMM,
  width = 12,
  height = 6,
  scale = 0.7
)

# Density plots - copied from Sofia's pipeline

log.cutoff <- log2(expr_cutoff)
nsamples <- ncol(quant_cDNA_DGE) 
col <- rainbow(nsamples) 
lcpm.Raw <- cpm(quant_cDNA_DGE$counts, log=TRUE) 
lcpm.Filt1 <- cpm(quant_cDNA_DGE_filter$counts, log=TRUE) 

png("./output/plots_QC/Density of count values - cDNA and ncRNA.png", width = 8, height = 4, units = 'in', res = 300) 
layout(matrix(1:2, nrow=1))

plot(density(lcpm.Raw[,1]),
     col=col[1], 
     xlim=c(-10,20),
     ylim=c(0,0.3),
     main="",
     xlab="") 

for (i in 2:nsamples){ 
  den <- density(lcpm.Raw[,i]) 
  lines(den$x, den$y, col=col[i])
} 

abline(v=log.cutoff,
       col="red",
       lwd=1,
       lty=2,
       main="") 

title("Raw data",
      xlab="log2-CPM")

plot(density(lcpm.Filt1[,1]),
     col=col[1],
     xlim=c(-10,20),
     ylim=c(0,0.3),
     main="",
     xlab="") 

for (i in 2:nsamples){
  den <- density(lcpm.Filt1[,i]) 
  lines(den$x, den$y, col=col[i])
} 

abline(v=log.cutoff,
       col="red",
       lwd=1,
       lty=2,
       main="") 

title("Filtered data (median CPM > 0.5)",
      xlab="log2-CPM") 

dev.off()

# Correlation heatmap - copied from Sofia's file

png(
  "./output/plots_QC/Sample correlation heatmap - cDNA.png",
  width = 10,
  height = 10,
  units = 'in',
  res = 300
)
par(mar = c(16, 16, 16, 16))

heatmap(cor(lcpm.Filt1))
title("Filtered data")

dev.off()


# Save files 

saveRDS(quant_cDNA_DGE_filter, file = "./output/data_expression/pre_DGE/quant_cDNA_DGE_filter.RDS")
