### Don't source this file by itself; call in from another file after running env_prep.R
### Load dataset 

quant_cDNA_DGE_clean <- readRDS(file = "./output/quant_cDNA_DGE_edgeRfilter.RDS")

quant_cDNA_ncRNA_ENSEMBL_DGE_clean <- readRDS(file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE_edgeRfilter.RDS")

### Make design matrix 

## Treat cell line as a random effect (duplicateCorrelation)
## Treat passage and treatment as fixed effects
## Treat batch effect as random intercept...?

table_design <- quant_cDNA_DGE_clean$samples %>%
  mutate(Passage = factor(Passage,
                             levels = c("P5", "P13"))) %>%
  mutate(Treatment = factor(Treatment, 
                               levels = c("Untreated", "Treated")))

design <- model.matrix(~run_date + Passage*Treatment,
                       data = table_design)

## Apply voom transformation

quant_cDNA_voom <- voom(quant_cDNA_DGE,
                            design,
                            plot = TRUE)

quant_cDNA_ncRNA_ENSEMBL_voom <- voom(quant_cDNA_ncRNA_ENSEMBL_DGE,
                            design,
                            plot = TRUE)

## MDS plot

plotMDS(quant_cDNA_voom,
        labels = table_design$run_date)


