here::i_am("R/7_post_GSVA_5_POI_boxplot.R")

########################
# Build boxplots for pathways of interest
########################

# Import packages
library(DESeq2)
library(GSVA)
library(here)
library(limma)
library(tidyverse)

# Load data
fit_gsva <- readRDS(
  here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "GSVA_results.RDS"
  )
)

quant_gsva <- readRDS(
  here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "quant_GSVA.RDS"
  )
)

# Set pathways and collections of interest
pathways_interest <- list(
  h = list(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_E2F_TARGETS"
  ),
  c2_cgp = list(
    "FOROUTAN_INTEGRATED_TGFB_EMT_UP",
    "AMIT_EGF_RESPONSE_40_HELA",
    "BROWNE_INTERFERON_RESPONSIVE_GENES",
    "VERNELL_RETINOBLASTOMA_PATHWAY_UP"
  ),
  c2_cp = list(
    "WP_CELL_DIFFERENTIATION_EXPANDED_INDEX",
    "WP_NUCLEAR_RECEPTORS",
    "WP_OVERVIEW_OF_INTERFERONS_MEDIATED_SIGNALING_PATHWAY",
    "REACTOME_DNA_METHYLATION"
  ),
  GOBP = list(
    "GOBP_CELL_FATE_DETERMINATION",
    "GOBP_INTEGRATED_STRESS_RESPONSE_SIGNALING",
    "GOBP_FAS_SIGNALING_PATHWAY",
    "GOBP_CELL_CYCLE_DNA_REPLICATION_INITIATION"
  )
)

fit_gsva_interest <- fit_gsva[c("h", "c2_cgp", "c2_cp", "GOBP")]
quant_gsva_interest <- quant_gsva[c("h", "c2_cgp", "c2_cp", "GOBP")]

# Define GSVA ES plotting functions
## For passage
plot_poi_passage <- function(fit, quant, geneset) {
  # Get table containing ES for gene set
  geneset_counts <- data.frame(
    score = exprs(quant)[geneset, ],
    Passage = phenoData(quant)$Passage,
    Treatment = phenoData(quant)$Treatment,
    Day = phenoData(quant)$Day
  ) %>%
    filter(
      Day == "D3",
      Treatment == "Untreated"
    )

  # Get position to start drawing signif bars
  y_position <- max(geneset_counts$score) * 1.1

  # Plot the data
  plot <- ggplot(
    geneset_counts,
    aes(
      x = Passage,
      y = score,
    )
  ) +
    geom_boxplot(
      aes(
        color = Passage,
        shape = Passage
      ),
      outliers = FALSE
    ) +
    geom_jitter(
      aes(
        fill = Passage,
        shape = Passage
      ),
      colour = "black",
      alpha = 0.7,
      stroke = 0.5
    ) +
    geom_signif(
      comparisons = list(
        c("P5", "P7"),
        c("P7", "P13"),
        c("P5", "P13")
      ),
      annotation = c(
        topTable(
          fit,
          coef = 13,
          sort.by = "none",
          number = Inf
        ) %>%
          .[geneset, "adj.P.Val"],
        topTable(
          fit,
          coef = 14,
          sort.by = "none",
          number = Inf
        ) %>%
          .[geneset, "adj.P.Val"],
        topTable(
          fit,
          coef = 15,
          sort.by = "none",
          number = Inf
        ) %>%
          .[geneset, "adj.P.Val"]
      ) %>%
        stars_pval(),
      y_position = c(y_position, y_position, y_position),
      textsize = 3,
      step_increase = 0.2
    ) +
    scale_colour_manual(
      values = palette_merge[c(1, 3, 5)]
    ) +
    scale_fill_manual(
      values = palette_merge[c(1, 3, 5)]
    ) +
    scale_shape_manual(
      values = c(21, 24, 22)
    ) +
    scale_y_continuous(expand = expansion(0, 0.3)) +
    theme_classic() +
    ggtitle(
      label = geneset %>%
        str_replace_all("HALLMARK", "H") %>%
        str_replace_all("_", " ") %>%
        str_wrap(
          width = 20,
          whitespace_only = TRUE
        )
    ) +
    ylab("GSVA ES")

  # Return object
  return(plot)
}

# Plot GSVA
## For passage
plots_poi_passage <- pmap(
  list(
    fit_gsva_interest,
    quant_gsva_interest,
    pathways_interest
  ),
  \(fit, quant, pathways_interest) {
    map(
      pathways_interest,
      \(geneset) {
        plot_poi_passage(fit, quant, geneset)
      }
    )
  }
)

# Save data
saveRDS(
  plots_poi_passage,
  here::here(
    "output",
    "data_enrichment",
    "GSVA",
    "plots_poi_passage.RDS"
  )
)
