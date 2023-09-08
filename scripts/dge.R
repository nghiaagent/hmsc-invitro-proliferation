### Don't source this file by itself; call in from another file after running env_prep.R
### Load dataset 
### Dataset has both mRNA and ncRNA

quant_DGE_clean <- readRDS(file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE_edgeRfilter.RDS")

### Make design matrix 

## Treat cell line as an additive factor (to use with duplicateCorrelation)
## Treat passage and treatment as fixed effects
## Treat batch effect as random intercept...?

table_design <- quant_DGE_clean$samples %>%
  mutate(Passage = factor(Passage,
                             levels = c("P5", "P13"))) %>%
  mutate(Treatment = factor(Treatment, 
                               levels = c("Untreated", "Treated")))

design <- model.matrix(~Passage*Treatment,
                       data = table_design)

## Apply voom transformation

quant_DGE_voom <- voom(quant_DGE_clean,
                            design,
                            plot = TRUE)

## MDS / PCA plot

Glimma::glMDSPlot(quant_DGE_clean,
                  labels = quant_DGE_clean$samples$Sample_ID,
                  groups = quant_DGE_clean$samples)



