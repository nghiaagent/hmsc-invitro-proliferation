# Looking closer into PC1 and PC2 of our data, which explains most of batch effect

source("./scripts/dge.R")


# Run PCA on top 90% variable genes

# Choose optimal number of PCs to retain
# Optimal number of PCs is between 1 - 4; very large PC1 due to batch fx

rownames(quant_DGE_clean$samples) <- quant_DGE_clean$samples$ID

pca <- PCAtools::pca(mat = quant_DGE_clean$counts,
                     metadata = quant_DGE_clean$samples,
                     removeVar = 0.1)

pca_optimisePCs <- PCAtools::parallelPCA(quant_DGE_clean$counts)

PCAtools::screeplot(pca,
          components = PCAtools::getComponents(pca, 1:30))


# Draw plot of PC1 and PC2, plus loadings, plus circles around batches

PCAtools::biplot(pca,
                 showLoadings = TRUE,
                 colby = "run_date",
                 encircle = TRUE,
                 encircleFill = TRUE,
                 title = 'Strong batch effect',
                 legendPosition = "bottom")

## Top PC1 contributors are ACTB, VIM, TMSB10. That's bad news.

# Draw paired biplot of first 6 PCs

## Batch effect - captured by PC1; PC2

PCAtools::pairsplot(pca,
                    components = PCAtools::getComponents(pca, 1:6),
                    triangle = TRUE,
                    colby = "run_date",
                    gridlines.minor = FALSE, 
                    hline = 0, vline = 0,
                    title = 'Pairs plot - grouped by batch',
                    plotaxes = FALSE)

## Difference between cell pops - captured by PC3

PCAtools::pairsplot(pca,
                    components = PCAtools::getComponents(pca, 1:6),
                    triangle = TRUE,
                    colby = "cell_line",
                    gridlines.minor = FALSE, 
                    hline = 0, vline = 0,
                    title = 'Pairs plot - grouped by cell population',
                    plotaxes = FALSE)

## Difference between passages - captured by PC5, 6; although very low variance.

PCAtools::pairsplot(pca,
                    components = PCAtools::getComponents(pca, 1:6),
                    triangle = TRUE,
                    colby = "Passage",
                    gridlines.minor = FALSE, 
                    hline = 0, vline = 0,
                    title = 'Pairs plot - grouped by growth phase (A vs. C)',
                    plotaxes = FALSE)

## Difference between treatments - not captured by these PCs

PCAtools::pairsplot(pca,
                    components = PCAtools::getComponents(pca, 1:6),
                    triangle = TRUE,
                    colby = "Treatment",
                    gridlines.minor = FALSE, 
                    title = 'Pairs plot - grouped by treatment',
                    plotaxes = FALSE)

## TL;DR:
## PC1, 2: Batch FX
## PC3: cell populations
## PC5, 6: passages
## No wonder it's hard to compare between treatments

# Determine genes driving variation in each PC

plot_loadings <- PCAtools::plotloadings(pca,
             rangeRetain = 0.0032,
             labSize = 4.0,
             title = 'Loadings plot',
             subtitle = 'Top 0.32% variables',
             shape = 24,
             col = c('limegreen', 'black', 'red3'),
             drawConnectors = TRUE)

## Top genes are ACTB, ACTG1, GAPDH, RMRP, H1RNA. They shouldn't be the most variable genes, right...

# Associate PCs to variables

plot_pcacor <- PCAtools::eigencorplot(pca,
                                      components = PCAtools::getComponents(pca, 1:10),
                                      metavars = c('Passage', 
                                                   'Treatment',
                                                   'condition_ID',
                                                   'cell_line',
                                                   'run_date'),
                                      col = c('darkblue', 'blue2', 'black', 'red2', 'darkred'),
                                      cexCorval = 0.7,
                                      colCorval = 'white',
                                      fontCorval = 2,
                                      posLab = 'bottomleft',
                                      rotLabX = 45,
                                      posColKey = 'top',
                                      cexLabColKey = 1.5,
                                      scale = TRUE,
                                      main = 'PC1-10 correlations with variables',
                                      colFrame = 'white',
                                      plotRsquared = FALSE)

## According to this, PC3 captures batch effect + cell pop differences; PC4 + 6 captures changes between growth phases. 
## None captures changes between treatments

# More biplots based on the correlation results

PCAtools::biplot(pca,
                 x = "PC3",
                 y = "PC1",
                 lab = NULL,
                 showLoadings = TRUE,
                 colby = "run_date",
                 encircle = TRUE,
                 encircleFill = TRUE,
                 title = 'Batch effect captured by PC 1 and 3',
                 legendPosition = "bottom")

PCAtools::biplot(pca,
                 x = "PC3",
                 y = "PC2",
                 lab = NULL,
                 showLoadings = TRUE,
                 colby = "cell_line",
                 encircle = TRUE,
                 encircleFill = TRUE,
                 title = 'Differences between cell populations captured by PC 2 and 3',
                 legendPosition = "bottom")

PCAtools::biplot(pca,
                 x = "PC6",
                 y = "PC4",
                 lab = NULL,
                 showLoadings = TRUE,
                 colby = "Passage",
                 encircle = TRUE,
                 encircleFill = TRUE,
                 title = 'Difference between growth phases captured by PC 4 and 6',
                 legendPosition = "bottom")

# Ex