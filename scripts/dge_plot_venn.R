# Run DGE first

## Overlap of DEGs between passages, at D3 UT

png(
  "./output/Venn D3 UT.png",
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

## Overlap of DEGs between passages, at D5 UT

png(
  "./output/Venn D5 UT.png",
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