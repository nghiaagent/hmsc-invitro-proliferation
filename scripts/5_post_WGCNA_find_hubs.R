# Load data

## GCN with all samples

load(file.path(
  ".",
  "output",
  "data_WGCNA",
  "WGCNA_allsamples",
  "GCN_output.Rdata"
))

# Derive module assignment

genes_sel <- labels2colors(GCN$net$colors) %in% c("turquoise", "red")
module_assign <- labels2colors(GCN$net$colors)

# Derive module eigengenes

module_eigengenes <- moduleEigengenes(GCN$E, labels2colors(GCN$net$colors))$eigengenes

# Screen for hub genes based on gene significance and module membership
# Within red and turquoise modules

## Calculate intramodular connectivity (degree) of each gene and whole network connectivity
### Intramodular connectivity: abs(cor(GCN$E, use = "pearson")) ^ 6
### kTotal: Whole network connectivity
### kWithin: Within module connectivity
### kOut: Outside module connectivity
### kDiff = kWithin - kOut
## Calculate gene significance of each gene in chosen modules

scores_connectivity <- intramodularConnectivity(abs(cor(GCN$E, use = "p")) ^
                                                  6, vec_modules) %>%
  cbind(as.numeric(cor(GCN$targets$Passage, GCN$E, use = "p"))) %>%
  rename("gene_significance" = 5) %>%
  cbind(module_assign)

## Calculate module membership for each gene

scores_membership <- signedKME(GCN$E, module_eigengenes, outputColumnName = "membership_")

## Combine all data

scores_hub <- cbind(scores_connectivity, scores_membership)

scores_hub_sel <- scores_hub[genes_sel, ]

## Plot GS against intramodular connectivity

plot_scores_connectivity <- ggscatter(
  data = scores_hub_sel,
  x = "kWithin",
  y = "gene_significance",
  color = "module_assign",
  add = "reg.line",
  conf.int = TRUE,
  alpha = 0.1
) +
  stat_cor() +
  scale_color_manual(values = c("red", "turquoise")) +
  labs(x = "Intramodular connectivity", y = "Gene significance for passage") +
  theme(legend.position = "none") +
  facet_wrap( ~ module_assign, labeller = labeller(
    module_assign = c("red" = "Module: red", "turquoise" = "Module: turquoise")
  ))

## Plot GS against module membership

plot_scores_membership_turquoise <- ggscatter(
  data = scores_hub_sel,
  x = "membership_turquoise",
  y = "gene_significance",
  color = "turquoise",
  add = "reg.line",
  conf.int = TRUE,
  alpha = 0.1
) +
  geom_vline(xintercept = 0.85) +
  geom_hline(yintercept = c(-0.4, 0.4)) +
  stat_cor() +
  labs(x = "Membership in turquoise module", y = "Gene significance for passage") +
  theme(legend.position = "none")

plot_scores_membership_red <- ggscatter(
  data = scores_hub_sel,
  x = "membership_red",
  y = "gene_significance",
  color = "red",
  add = "reg.line",
  conf.int = TRUE,
  alpha = 0.1
) +
  geom_vline(xintercept = 0.85) +
  geom_hline(yintercept = c(-0.4, 0.4)) +
  stat_cor() +
  labs(x = "Membership in red module", y = "Gene significance for passage") +
  theme(legend.position = "none")

grid <- plot_grid(
  plot_scores_connectivity,
  plot_grid(
    plot_scores_membership_red,
    plot_scores_membership_turquoise,
    ncol = 2,
    nrow = 1,
    labels = c("B", "C")
  ),
  ncol = 1,
  nrow = 2,
  labels = c("A", "")
)

## Screen based on GS and MM metrics in each module

scores_hub_red <- filter(scores_hub,
                         abs(gene_significance) > 0.4,
                         membership_red > 0.85,
                         module_assign == "red") %>%
  select(gene_significance, membership_red, kWithin) %>%
  mutate(ENTREZID = rownames(.) %>% as.integer()) %>%
  left_join(GCN$genes)

scores_hub_turquoise <- filter(
  scores_hub,
  abs(gene_significance) > 0.4,
  membership_turquoise > 0.85,
  module_assign == "turquoise"
) %>%
  select(gene_significance, membership_turquoise, kWithin) %>%
  mutate(ENTREZID = rownames(.) %>% as.integer()) %>%
  left_join(GCN$genes)

# Export data

write_xlsx(
  list(
    "Turquoise hub candidates" = scores_hub_turquoise,
    "Red hub candidates" = scores_hub_red
  ),
  path = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "hub_candidates.xlsx"
  )
)

ggsave(
  filename = file.path(
    ".",
    "output",
    "plots_WGCNA",
    "WGCNA_allsamples",
    "Hub gene screening.png"
  ),
  grid,
  width = 9,
  height = 9
)

fwrite(
  scores_hub_turquoise$GENENAME %>% list(),
  file = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "turquoise_hub_candidates.txt"
  ),
  na = ""
)

fwrite(
  scores_hub_red$GENENAME %>% list(),
  file = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "red_hub_candidates.txt"
  ),
  na = ""
)