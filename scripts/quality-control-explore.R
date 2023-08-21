### Don't source this file by itself; call in from another file after running env_prep.R
### Load dataset 

quant_cDNA_DGE <- readRDS(file = "./output/quant_cDNA_DGE.RDS")

quant_cDNA_ncRNA_ENSEMBL_DGE <- readRDS(file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE.RDS")

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
       title = "Size of libraries mapped to different transcriptomes",
       fill = "Transcriptome type")

ggsave("./output/plots_QC/Library size.png",
       plot = plot_libsize)

## Crazy: Paired t-test to see if library size are 
## significantly different between mapped transcriptomes (cDNA only vs cDNA-ncRNA)

t.test(quant_libsize$lib.size_cDNA,
       quant_libsize$lib.size_cDNA_ncRNA_ENSEMBL,
       paired = TRUE)

## They are significantly different. More counts (95% CI 201143-304042) in the cDNA-ncRNA transcriptome. 

# Get average count for each gene in all samples

df_stats_cDNA <- rowMeans(quant_cDNA_DGE$counts) %>%
  as_tibble() %>%
  mutate(ENSEMBL = rownames(quant_cDNA_DGE$counts)) %>%
  relocate(ENSEMBL) %>%
  dplyr::rename("average_counts" = "value") %>%
  mutate(ratio_total_counts = average_counts/sum(average_counts)) %>%
  left_join(quant_cDNA_DGE$genes, by = join_by(ENSEMBL == ENSEMBL))


df_avg_counts <- tibble(quant_cDNA_DGE$genes$SYMBOL,
                        quant_cDNA_DGE$genes$GENENAME)