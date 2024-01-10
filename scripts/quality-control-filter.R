### Don't source this file by itself; call in from another file after running env_prep.R
### Load dataset 

quant_cDNA_ncRNA_ENSEMBL_DGE <- readRDS(file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE.RDS")

# Library size barplot. Library size is the total number of counts across all genes per sample

## Generate table for this purpose

quant_libsize <- dplyr::select(quant_cDNA_ncRNA_ENSEMBL_DGE$samples,
                                    c("lib.size", "ID")) %>%
  relocate(ID)

summary(quant_libsize$lib.size)

## Barplot

plot_libsize <- ggplot(quant_libsize,
       aes(x = ID,
           y = lib.size/1e6)
       ) +
  geom_bar(stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "",
       y = "Library size (million)",
       title = "Size of libraries")

ggsave("./output/plots_QC/Library size.png",
       plot = plot_libsize)

# Get average count for each gene in all samples, sorted from highest to lowest

df_stats_cDNA_ncRNA_ENSEMBL <- rowMeans(quant_cDNA_ncRNA_ENSEMBL_DGE$counts) %>%
  as_tibble() %>%
  mutate(ENSEMBL = rownames(quant_cDNA_ncRNA_ENSEMBL_DGE$counts)) %>%
  relocate(ENSEMBL) %>%
  dplyr::rename("average_counts" = "value") %>%
  mutate(ratio_total_counts = average_counts/sum(average_counts)) %>%
  left_join(quant_cDNA_ncRNA_ENSEMBL_DGE$genes, by = join_by(ENSEMBL == GENEID)) %>%
  arrange(desc(ratio_total_counts))

## Write list of top 50 genes

head(df_stats_cDNA_ncRNA_ENSEMBL, 50) %>%
  write_csv("./output/QC/Top_20_Expressed_Genes_cDNA_ncRNA_ENSEMBL.csv")

## Filter genes with low expression using edgeR function
# 
# expr_cutoff <- 40
# 
# keep_cDNA_ncRNA_ENSEMBL_edgeRfiter <- filterByExpr(quant_cDNA_ncRNA_ENSEMBL_DGE,
#                                      group = "condition_ID",
#                                      min.count = expr_cutoff,
#                                      min.total.count = 1.5*expr_cutoff)
# 
# quant_cDNA_ncRNA_ENSEMBL_DGE_edgeRfilter <- quant_cDNA_ncRNA_ENSEMBL_DGE[keep_cDNA_ncRNA_ENSEMBL_edgeRfiter, , keep.lib.sizes=FALSE] %>%
#   calcNormFactors()

## Filter genes based on CPM 0.5 threshold

expr_cutoff <- 0.5 # in cpm sum(median_cpm > expr_cutoff)

median_cpm <- apply(cpm(quant_cDNA_ncRNA_ENSEMBL_DGE), 1, median)
quant_cDNA_ncRNA_ENSEMBL_DGE_filter <- quant_cDNA_ncRNA_ENSEMBL_DGE[median_cpm > expr_cutoff, ] %>%
  calcNormFactors()

# Per sample distribution; before and after adding TMM scaling factor

png("./output/plots_QC/Per sample counts distribution.png", width = 40, height = 60, units = 'cm', res = 600) 
par(mfrow=c(2,1))

lcpm <- cpm(quant_cDNA_ncRNA_ENSEMBL_DGE[apply(cpm(quant_cDNA_ncRNA_ENSEMBL_DGE), 1, median) > expr_cutoff, ],
            log=TRUE)

boxplot(lcpm, 
        las=2, 
        main="")

title(main="Unnormalised data",ylab="Log-cpm")

lcpm <- cpm(quant_cDNA_ncRNA_ENSEMBL_DGE_filter,
            log=TRUE)

boxplot(lcpm, 
        las=2,
        main="")

title(main="TMM-Normalised data",ylab="Log-cpm")

dev.off()

# Density plots - copied from Sofia's pipeline

log.cutoff <- log2(expr_cutoff)

png("./output/plots_QC/Density of count values - cDNA and ncRNA.png", width = 10, height = 30, units = 'cm', res = 600) 
nsamples <- ncol(quant_cDNA_ncRNA_ENSEMBL_DGE) 
col <- rainbow(nsamples) 
par(mfrow=c(2,1)) 
lcpm.Raw <- cpm(quant_cDNA_ncRNA_ENSEMBL_DGE$counts, log=TRUE) 
plot(density(lcpm.Raw[,1]), col=col[1], 
     xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="") 
for (i in 2:nsamples){ den <-density(lcpm.Raw[,i]) 
lines(den$x, den$y, col=col[i]) } 
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="") 
title("Raw data",xlab="log2-CPM")

lcpm.Filt1 <- cpm(quant_cDNA_ncRNA_ENSEMBL_DGE_filter$counts, log=TRUE) 
plot(density(lcpm.Filt1[,1]), col=col[1], xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="") 
for (i in 2:nsamples){ den <-density(lcpm.Filt1[,i]) 
lines(den$x, den$y, col=col[i]) } 
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="") 
title("Filtered data (median CPM > 0.5)",xlab="log2-CPM") 

dev.off()

# Heatmap - copied from Sofia's file

png(
  "./output/plots_QC/Sample correlation heatmap - cDNA and ncRNA.png",
  width = 60,
  height = 60,
  units = 'cm',
  res = 300
)
par(mfrow = c(1, 2))
lcpm.Raw <- cpm(quant_cDNA_ncRNA_ENSEMBL_DGE$counts, log = TRUE)
heatmap(cor(lcpm.Raw))
title("Raw data")
lcpm.Filt <- cpm(quant_cDNA_ncRNA_ENSEMBL_DGE_filter$counts, log = TRUE)
heatmap(cor(lcpm.Filt))
title("Filtered data")
dev.off()


# Save files 

saveRDS(quant_cDNA_ncRNA_ENSEMBL_DGE_filter, file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE_filter.RDS")
