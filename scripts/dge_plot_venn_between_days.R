# Load data

source("./scripts/dge_cellpops_as_fixed.R")

png(
  "./output/Venn T vs UT @ Passage.png",
  width = 15,
  height = 30,
  units = 'cm',
  res = 400
)

par(mfrow = c(3, 1))

vennDiagram(decideTests(fit_contrasts)[,c(7,10)],
            include = c("up",
                        "down"),
            names = c("Untreated",
                      "Treated"),
            counts.col = c("red",
                           "blue"))



vennDiagram(decideTests(fit_contrasts)[,c(8,11)],
            include = c("up",
                        "down"),
            names = c("Untreated",
                      "Treated"),
            counts.col = c("red",
                           "blue"))

vennDiagram(decideTests(fit_contrasts)[,c(9,12)],
                       include = c("up",
                                   "down"),
                       names = c("Untreated",
                                 "Treated"),
                       counts.col = c("red",
                                      "blue"))

dev.off()


png(
  "./output/Venn T vs UT @ D3.png",
  width = 15,
  height = 30,
  units = 'cm',
  res = 400
)

par(mfrow = c(3, 1))

vennDiagram(decideTests(fit_contrasts)[,c(13,16)],
            include = c("up",
                        "down"),
            names = c("Untreated",
                      "Treated"),
            counts.col = c("red",
                           "blue"))



vennDiagram(decideTests(fit_contrasts)[,c(14,17)],
            include = c("up",
                        "down"),
            names = c("Untreated",
                      "Treated"),
            counts.col = c("red",
                           "blue"))

vennDiagram(decideTests(fit_contrasts)[,c(15,18)],
            include = c("up",
                        "down"),
            names = c("Untreated",
                      "Treated"),
            counts.col = c("red",
                           "blue"))

dev.off()


png(
  "./output/Venn T vs UT @ D5.png",
  width = 15,
  height = 30,
  units = 'cm',
  res = 400
)

par(mfrow = c(3, 1))

vennDiagram(decideTests(fit_contrasts)[,c(19,22)],
            include = c("up",
                        "down"),
            names = c("Untreated",
                      "Treated"),
            counts.col = c("red",
                           "blue"))



vennDiagram(decideTests(fit_contrasts)[,c(20,23)],
            include = c("up",
                        "down"),
            names = c("Untreated",
                      "Treated"),
            counts.col = c("red",
                           "blue"))

vennDiagram(decideTests(fit_contrasts)[,c(21,24)],
            include = c("up",
                        "down"),
            names = c("Untreated",
                      "Treated"),
            counts.col = c("red",
                           "blue"))

dev.off()