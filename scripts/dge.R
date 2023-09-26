### Don't source this file by itself; call in from another file after running env_prep.R
# Load dataset 
## Dataset has both mRNA and ncRNA

quant_DGE_clean <- readRDS(file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE_edgeRfilter.RDS")

## MDS / PCA plot

Glimma::glMDSPlot(quant_DGE_clean,
                  labels = quant_DGE_clean$samples$Sample_ID,
                  groups = quant_DGE_clean$samples)

### Make design matrix 

## Treat cell line as an additive factor
## Treat passage and treatment as fixed effects

table_design <- quant_DGE_clean$samples %>%
  mutate(Passage = factor(Passage,
                             levels = c("P5", "P7", "P13"))) %>%
  mutate(Treatment = factor(Treatment, 
                               levels = c("Untreated", "Treated")))

design <- model.matrix(~Passage*Treatment + Passage*Treatment:cell_line,
                       data = table_design)

# Use ComBat-seq to correct for batch fx
# Covariates considered: Passage, Treatment, cell population

quant_DGE_batchcor_withcovariates <-
  sva::ComBat_seq(quant_DGE_clean$counts,
                  batch = table_design$run_date,
                  covar_mod = design)

quant_DGE_clean_batchcor <- quant_DGE_clean
quant_DGE_clean_batchcor$counts <- quant_DGE_batchcor_withcovariates

# Subset data to just P5 vs P13; D3

quant_DGE_clean_batchcor_subset <- quant_DGE_clean_batchcor[,which(
  quant_DGE_clean_batchcor$samples$Day == "D3" &
    !quant_DGE_clean_batchcor$samples$Passage == "P7"
)]

table_design <- quant_DGE_clean_batchcor_subset$samples %>%
  mutate(Passage = factor(Passage,
                          levels = c("P5", "P13"))) %>%
  mutate(Treatment = factor(Treatment, 
                            levels = c("Untreated", "Treated")))

design <- model.matrix(~Passage*Treatment + Passage*Treatment:cell_line,
                       data = table_design)

## Apply voom transformation
### And output mean-variance trend plot

png("./output/plots_QC/voom mean-variance trend.png", width = 20, height = 20, units = 'cm', res = 600) 

quant_DGE_voom <- voom(quant_DGE_clean_batchcor_subset,
                            design,
                            plot = TRUE)

dev.off()


## Duplicate correlation not necessary because cell_line is treated as a fixed (additive) effect

fit <- lmFit(quant_DGE_voom,
             design) %>%
  eBayes()

summary(decideTests(fit))

## Draw interactive Volcano plot
## coef chooses column of design matrix (factor) that the plot shows

Glimma::glimmaVolcano(fit,
                      quant_DGE_voom,
                      coef = 3)

# Define posthoc tests
## Between treatments, at each passage, at each cell population
