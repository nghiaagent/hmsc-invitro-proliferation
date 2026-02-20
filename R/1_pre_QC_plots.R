# Load data
quant_cDNA_deseq <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "pre_DGE",
    "quant_cDNA_deseq.RDS"
  )
)
quant_cDNA_deseq_nofilter <- readRDS(
  file = here::here(
    "output",
    "data_expression",
    "pre_DGE",
    "quant_cDNA_deseq_nofilter.RDS"
  )
)

coldata <- colData(quant_cDNA_deseq) %>%
  as.data.frame()

x <- as.DGEList(quant_cDNA_deseq) %>%
  normLibSizes() %>%
  cpm(normalized.lib.sizes = FALSE)

# Get data for plotting libsize
list_counts <- list(
  norm = quant_cDNA_deseq %>%
    estimateSizeFactors() %>%
    counts(normalized = TRUE),
  nonorm = quant_cDNA_deseq %>%
    counts(normalized = FALSE),
  norm_TMM = as.DGEList(quant_cDNA_deseq) %>%
    normLibSizes() %>%
    cpm(normalized.lib.sizes = FALSE),
  nofilter_norm = quant_cDNA_deseq_nofilter %>%
    estimateSizeFactors() %>%
    counts(normalized = TRUE),
  nofilter_nonorm = quant_cDNA_deseq_nofilter %>%
    counts(normalized = FALSE)
) %>%
  imap(
    \(counts, name) {
      counts %>%
        as.data.frame() %>%
        rownames_to_column(var = "gene") %>%
        pivot_longer(
          cols = !gene,
          names_to = "sample_id",
          values_to = "counts"
        ) %>%
        mutate(logcounts_offset = log10(counts + 1)) %>%
        mutate(data_type = name) %>%
        left_join(
          coldata,
          by = join_by("sample_id" == "names")
        )
    }
  ) %>%
  rbindlist()

# Calculate total library size
df_libsize <- list_counts %>%
  filter(data_type == "nofilter_nonorm") %>%
  group_by(sample_id) %>%
  summarise(libsize = sum(counts)) %>%
  left_join(
    coldata,
    by = join_by("sample_id" == "names")
  ) %>%
  arrange(cell_line) %>%
  arrange(
    Passage,
    Day,
    Treatment
  ) %>%
  mutate(
    sample_id = sample_id %>%
      factor() %>%
      fct_inorder()
  )

# Plot libsize
plot_libsize <- df_libsize %>%
  ggplot(
    aes(
      x = sample_id,
      y = libsize / 1000000,
      fill = timepoint_ID
    )
  ) +
  geom_col(
    colour = "white"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title.x = element_blank(),
    legend.title = element_text()
  ) +
  scale_fill_manual(
    values = palette_merge[1:6]
  ) +
  geom_hline(
    yintercept = c(
      5,
      10,
      15
    ),
    linetype = "dashed",
    color = "grey30"
  ) +
  labs(
    y = "Library size (million)",
    fill = "Timepoint"
  )

# Plot boxplots of normalised vs. unnormalised data
plots_normalisation <- list_counts %>%
  filter(data_type %in% c("nonorm", "norm", "norm_TMM")) %>%
  mutate(
    sample_id = sample_id %>%
      factor(levels = levels(df_libsize$sample_id))
  ) %>%
  ggplot(
    aes(
      x = sample_id,
      y = logcounts_offset,
      colour = timepoint_ID
    )
  ) +
  geom_boxplot() +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    axis.title.x = element_blank(),
    legend.title = element_text(),
    legend.position = "none"
  ) +
  scale_colour_manual(
    values = palette_merge[1:6]
  ) +
  labs(
    y = "log10(counts)"
  ) +
  facet_wrap(
    ~data_type,
    labeller = as_labeller(
      c(
        "nonorm" = "Unnormalised data",
        "norm" = "DESeq2 normalised data",
        "norm_TMM" = "TMM normalised data"
      )
    ),
    ncol = 1,
    scales = "free_y"
  )

# Plot density of genes before and after normalisation
plots_filtering <- list_counts %>%
  filter(data_type %in% c("nofilter_norm", "norm")) %>%
  ggplot(
    aes(
      x = logcounts_offset,
      colour = timepoint_ID
    )
  ) +
  geom_density(
    linewidth = 0.5,
    alpha = 0.6
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    ),
    legend.title = element_text(),
    legend.position = "none"
  ) +
  scale_colour_manual(
    values = palette_merge[1:6]
  ) +
  geom_vline(
    xintercept = log10(41),
    linetype = "dashed",
    colour = "grey30"
  ) +
  labs(
    x = "log10(counts)",
    y = "Density"
  ) +
  facet_wrap(
    ~data_type,
    ncol = 2,
    labeller = as_labeller(
      c(
        "nofilter_norm" = "Unfiltered data",
        "norm" = "FilterByExpr 40 counts"
      )
    ),
    scales = "free_y"
  )

# Save plots
## Library size
ggsave(
  here::here(
    "output",
    "plots_QC",
    "Library size.png"
  ),
  plot_libsize,
  width = 16,
  height = 9,
  scale = 0.7
)

## Distribution
ggsave(
  here::here(
    "output",
    "plots_QC",
    "Distribution.png"
  ),
  plots_normalisation,
  width = 8,
  height = 10,
  scale = 1
)

## Gene filtering
ggsave(
  here::here(
    "output",
    "plots_QC",
    "Gene filter.png"
  ),
  plots_filtering,
  width = 8,
  height = 6,
  scale = 0.7
)
