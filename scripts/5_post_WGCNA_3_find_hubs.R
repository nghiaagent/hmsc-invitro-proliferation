# Analyse modules for highly connected genes
## Turquoise: Positive correlation, day + passage
## Pink:      Positive correlation, day only
## Blue:      Negative correlation, day + passage
## Red:       Negative correlation, day only

# Load data
gcn <- readRDS(
  file = here::here(
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "GCN_output.RDS"
  )
)

# Set modules
modules_sel <- c(
  "turquoise",
  "pink",
  "blue",
  "red"
)

# Set thresholds

threshold_mm <- 0.8
threshold_gs <- 0.4

# Derive module assignment
genes_sel <- labels2colors(gcn$net$colors) %in% modules_sel
module_assign <- labels2colors(gcn$net$colors)

# Derive module eigengenes
module_eigengenes <- moduleEigengenes(
  gcn$E,
  labels2colors(gcn$net$colors)
) %>%
  .$eigengenes

# Screen for hub genes based on gene significance and module membership
# Within selected modules
# Associate with passage changes

## Calculate intramodular connectivity (degree) of each gene and whole network connectivity
### Intramodular connectivity: abs(cor(GCN$E, use = "pearson")) ^ 6
### kTotal: Whole network connectivity
### kWithin: Within module connectivity
### kOut: Outside module connectivity
### kDiff = kWithin - kOut
## Calculate gene significance of each gene in chosen modules

scores_connectivity <- intramodularConnectivity(
  abs(cor(gcn$E, use = "p"))^6,
  module_assign
) %>%
  cbind(as.numeric(
    cor(
      gcn$targets$Passage,
      gcn$E,
      use = "p"
    )
  )) %>%
  rename("gene_significance" = 5) %>%
  cbind(module_assign)

## Calculate module membership for each gene
scores_membership <- signedKME(
  gcn$E,
  module_eigengenes,
  outputColumnName = "membership_"
)

## Combine all data
scores_hub <- cbind(
  scores_connectivity,
  scores_membership
)

scores_hub_sel <- scores_hub[genes_sel, ] %>%
  mutate(
    module_assign = factor(
      module_assign,
      levels = modules_sel
    )
  )

## Plot GS against intramodular connectivity
plot_scores_connectivity <- ggscatter(
  data = scores_hub_sel,
  x = "kWithin",
  y = "gene_significance",
  color = "module_assign",
  add = "reg.line",
  conf.int = TRUE,
  alpha = 0.1,
  size = 1
) +
  stat_cor() +
  scale_color_manual(values = modules_sel) +
  labs(
    x = "Intramodular connectivity",
    y = "Gene significance for passage"
  ) +
  theme(legend.position = "none") +
  facet_wrap(
    ~module_assign,
    ncol = 4,
    nrow = 1
  )

## Plot GS against module membership

plots_module_membership <- imap(
  list(
    "turquoise" = "membership_turquoise",
    "pink" = "membership_pink",
    "blue" = "membership_blue",
    "red" = "membership_red"
  ),
  \(x, color) {
    ggscatter(
      data = scores_hub_sel,
      x = x,
      y = "gene_significance",
      color = color,
      add = "reg.line",
      conf.int = TRUE,
      alpha = 0.1,
      size = 1
    ) +
      geom_vline(xintercept = threshold_mm) +
      geom_hline(yintercept = c(-1 * threshold_gs, threshold_gs)) +
      stat_cor() +
      labs(
        x = str_c(
          "Membership in ",
          color,
          " module"
        ),
        y = "Gene significance for passage"
      ) +
      theme(legend.position = "none")
  }
)

grid <- plot_grid(
  plot_scores_connectivity,
  plot_grid(
    plotlist = plots_module_membership,
    ncol = 4,
    nrow = 1,
    labels = c(
      "B",
      "",
      "",
      ""
    )
  ),
  ncol = 1,
  nrow = 2,
  labels = c(
    "A",
    ""
  )
)

## Screen based on GS and MM metrics in each module
scores_hub_turquoise <- filter(
  scores_hub,
  abs(gene_significance) > threshold_gs,
  membership_turquoise > threshold_mm,
  module_assign == "turquoise"
) %>%
  select(gene_significance, membership_turquoise, kWithin) %>%
  mutate(entrezid = rownames(.)) %>%
  left_join(
    rownames_to_column(gcn$genes),
    by = join_by("entrezid" == "rowname")
  )

scores_hub_pink <- filter(
  scores_hub,
  abs(gene_significance) > threshold_gs,
  membership_pink > threshold_mm,
  module_assign == "pink"
) %>%
  select(gene_significance, membership_pink, kWithin) %>%
  mutate(entrezid = rownames(.)) %>%
  left_join(
    rownames_to_column(gcn$genes),
    by = join_by("entrezid" == "rowname")
  )

scores_hub_blue <- filter(
  scores_hub,
  abs(gene_significance) > threshold_gs,
  membership_blue > threshold_mm,
  module_assign == "blue"
) %>%
  select(gene_significance, membership_blue, kWithin) %>%
  mutate(entrezid = rownames(.)) %>%
  left_join(
    rownames_to_column(gcn$genes),
    by = join_by("entrezid" == "rowname")
  )

scores_hub_red <- filter(
  scores_hub,
  abs(gene_significance) > threshold_gs,
  membership_red > threshold_mm,
  module_assign == "red"
) %>%
  select(gene_significance, membership_red, kWithin) %>%
  mutate(entrezid = rownames(.)) %>%
  left_join(
    rownames_to_column(gcn$genes),
    by = join_by("entrezid" == "rowname")
  )

# Export data
write.xlsx(
  list(
    "Turquoise hub candidates" = scores_hub_turquoise,
    "Red hub candidates" = scores_hub_red,
    "Blue hub candidates" = scores_hub_blue,
    "Pink hub candidates" = scores_hub_pink
  ),
  file = here::here(
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "hub_candidates.xlsx"
  )
)

ggsave(
  filename = here::here(
    "output",
    "plots_WGCNA",
    "WGCNA_allsamples",
    "Hub gene screening.png"
  ),
  grid,
  width = 12,
  height = 6
)

imap(
  list(
    "Turquoise hub candidates" = scores_hub_turquoise,
    "Pink hub candidates" = scores_hub_pink,
    "Blue hub candidates" = scores_hub_blue,
    "Red hub candidates" = scores_hub_red
  ),
  \(table, name) {
    fwrite(
      list(table$symbol),
      file = here::here(
        "output",
        "data_WGCNA",
        "WGCNA_allsamples",
        str_c(name, ".txt")
      ),
      na = ""
    )
  }
)
