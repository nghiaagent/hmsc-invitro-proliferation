### Don't source this file by itself; call in from another file after running env_prep.R
### Load dataset 

quant_cDNA_DGE <- readRDS(file = "./output/quant_cDNA_DGE.RDS")

quant_cDNA_ncRNA_ENSEMBL_DGE <- readRDS(file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE.RDS")

# Remove samples that are extreme outliers

list_remove <- c("hMSC-20176_P13_D5_0-3",
                 "hMSC-21558_P5_D3_0-2",
                 "hMSC-21558_P5_D3_10-3")

quant_cDNA_DGE <- quant_cDNA_DGE[,which(
  !quant_cDNA_DGE$samples$ID %in% list_remove
)]

quant_cDNA_ncRNA_ENSEMBL_DGE <- quant_cDNA_ncRNA_ENSEMBL_DGE[,which(
  !quant_cDNA_ncRNA_ENSEMBL_DGE$samples$ID %in% list_remove
)]

# Library size barplot. Library size is the total number of counts across all genes per sample

## Generate table for this purpose

quant_cDNA_libsize <- dplyr::select(quant_cDNA_DGE$samples,
                                    c("lib.size", "ID")) %>%
  dplyr::rename("lib.size_cDNA" = "lib.size")

quant_cDNA_ncRNA_ENSEMBL_libsize <- dplyr::select(quant_cDNA_ncRNA_ENSEMBL_DGE$samples,
                                    c("lib.size", "ID")) %>%
  dplyr::rename("lib.size_cDNA_ncRNA_ENSEMBL" = "lib.size")

quant_libsize <- full_join(quant_cDNA_libsize, quant_cDNA_ncRNA_ENSEMBL_libsize) %>%
  relocate(ID)

quant_libsize_plot <- pivot_longer(quant_libsize,
                                   !ID,
                                   names_to = "transcriptome_type",
                                   names_prefix = "lib.size_",
                                   values_to = "lib.size")

summary(quant_libsize$lib.size_cDNA)

summary(quant_libsize$lib.size_cDNA_ncRNA_ENSEMBL)

## Barplot

plot_libsize <- ggplot(quant_libsize_plot,
       aes(x = ID,
           y = lib.size/1e6,
           fill = transcriptome_type)
       ) +
  geom_bar(position = "dodge",
           stat = "identity") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = "",
       y = "Library size (million)",
       title = "Size of libraries mapped to mRNA only vs mRNA + ncRNA transcriptomes",
       fill = "Transcriptome type")

ggsave("./output/plots_QC/Library size.png",
       plot = plot_libsize)

## Crazy: Paired t-test to see if library size are 
## significantly different between mapped transcriptomes (cDNA only vs cDNA-ncRNA)

t.test(quant_libsize$lib.size_cDNA,
       quant_libsize$lib.size_cDNA_ncRNA_ENSEMBL,
       paired = TRUE)

## They are significantly different. More counts (95% CI 201143-304042) in the cDNA-ncRNA transcriptome. 

# Get average count for each gene in all samples, sorted from highest to lowest

df_stats_cDNA <- rowMeans(quant_cDNA_DGE$counts) %>%
  as_tibble() %>%
  mutate(ENSEMBL = rownames(quant_cDNA_DGE$counts)) %>%
  relocate(ENSEMBL) %>%
  dplyr::rename("average_counts" = "value") %>%
  mutate(ratio_total_counts = average_counts/sum(average_counts)) %>%
  left_join(quant_cDNA_DGE$genes, by = join_by(ENSEMBL == GENEID)) %>%
  arrange(desc(ratio_total_counts))

df_stats_cDNA_ncRNA_ENSEMBL <- rowMeans(quant_cDNA_ncRNA_ENSEMBL_DGE$counts) %>%
  as_tibble() %>%
  mutate(ENSEMBL = rownames(quant_cDNA_ncRNA_ENSEMBL_DGE$counts)) %>%
  relocate(ENSEMBL) %>%
  dplyr::rename("average_counts" = "value") %>%
  mutate(ratio_total_counts = average_counts/sum(average_counts)) %>%
  left_join(quant_cDNA_ncRNA_ENSEMBL_DGE$genes, by = join_by(ENSEMBL == GENEID)) %>%
  arrange(desc(ratio_total_counts))

## Write list of top 50 genes

head(df_stats_cDNA, 50) %>%
  write_csv("./output/QC/Top_20_Expressed_Genes_cDNA.csv")


head(df_stats_cDNA_ncRNA_ENSEMBL, 50) %>%
  write_csv("./output/QC/Top_20_Expressed_Genes_cDNA_ncRNA_ENSEMBL.csv")

## Filter genes with low expression using edgeR function

expr_cutoff <- 40

keep_cDNA_edgeRfiter <- filterByExpr(quant_cDNA_DGE,
                                      group = "condition_ID",
                                     min.count = expr_cutoff,
                                     min.total.count = 60)

quant_cDNA_DGE_edgeRfilter <- quant_cDNA_DGE[keep_cDNA_edgeRfiter, , keep.lib.sizes=FALSE] %>%
  calcNormFactors()

keep_cDNA_ncRNA_ENSEMBL_edgeRfiter <- filterByExpr(quant_cDNA_ncRNA_ENSEMBL_DGE,
                                     group = "condition_ID",
                                     min.count = expr_cutoff,
                                     min.total.count = 60)

quant_cDNA_ncRNA_ENSEMBL_DGE_edgeRfilter <- quant_cDNA_ncRNA_ENSEMBL_DGE[keep_cDNA_ncRNA_ENSEMBL_edgeRfiter, , keep.lib.sizes=FALSE] %>%
  calcNormFactors()

## cDNA only: from 38254 genes to 13780 genes (out of 18k genes targeted by AmpliSeq)
## cDNA + ncRNA from ENSEMBL: from 64541 to 14457 genes (out of 18k + 2k genes targeted by AmpliSeq)

# Heatmap and density plots - copied from Sofia's pipeline

log.cutoff <- log2(expr_cutoff)

## cDNA only

png("./output/plots_QC/Density of count values - cDNA only.png", width = 10, height = 30, units = 'cm', res = 600) 
nsamples <- ncol(quant_cDNA_DGE) 
col <- rainbow(nsamples) 
par(mfrow=c(2,1)) 
lcpm.Raw <- cpm(quant_cDNA_DGE$counts, log=TRUE) 
plot(density(lcpm.Raw[,1]), col=col[1], 
     xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="") 
for (i in 2:nsamples){ den <-density(lcpm.Raw[,i]) 
lines(den$x, den$y, col=col[i]) } 
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="") 
title("Raw data",xlab="log2-CPM")

lcpm.Filt1 <- cpm(quant_cDNA_DGE_edgeRfilter$counts, log=TRUE) 
plot(density(lcpm.Filt1[,1]), col=col[1], xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="") 
for (i in 2:nsamples){ den <-density(lcpm.Filt1[,i]) 
lines(den$x, den$y, col=col[i]) } 
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="") 
title("Filtered data (edgeR FilterByExpr())",xlab="log2-CPM") 

dev.off()

## cDNA and ncRNA

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

lcpm.Filt1 <- cpm(quant_cDNA_ncRNA_ENSEMBL_DGE_edgeRfilter$counts, log=TRUE) 
plot(density(lcpm.Filt1[,1]), col=col[1], xlim=c(-10,20), ylim=c(0,0.3), main="", xlab="") 
for (i in 2:nsamples){ den <-density(lcpm.Filt1[,i]) 
lines(den$x, den$y, col=col[i]) } 
abline(v=log.cutoff, col="red", lwd=1, lty=2, main="") 
title("Filtered data (edgeR FilterByExpr())",xlab="log2-CPM") 

dev.off()

# Save files 

saveRDS(quant_cDNA_DGE_edgeRfilter, file = "./output/quant_cDNA_DGE_edgeRfilter.RDS")

saveRDS(quant_cDNA_ncRNA_ENSEMBL_DGE_edgeRfilter, file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE_edgeRfilter.RDS")
