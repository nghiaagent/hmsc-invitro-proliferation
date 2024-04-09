# Load data

source("./scripts/dge_cellpops_as_fixed.R")

png(
  "./output/Venn D3.png",
  width = 30,
  height = 20,
  units = 'cm',
  res = 400
)

vennDiagram(decideTests(fit_contrasts)[,13:15],
            include = c("up",
                        "down"),
            names = c("Phase B - Phase A",
                      "Phase C - Phase B",
                      "Phase C - Phase A"),
            counts.col = c("red",
                           "blue"))

dev.off()


png(
  "./output/Venn D5.png",
  width = 30,
  height = 20,
  units = 'cm',
  res = 400
)

vennDiagram(decideTests(fit_contrasts)[,19:21],
            include = c("up",
                        "down"),
            names = c("Phase B - Phase A",
                      "Phase C - Phase B",
                      "Phase C - Phase A"),
            counts.col = c("red",
                           "blue"))

dev.off()