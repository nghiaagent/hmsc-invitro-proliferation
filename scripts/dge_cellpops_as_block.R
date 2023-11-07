### Don't source this file by itself; call in from another file after running env_prep.R
# Load dataset
## Dataset has both mRNA and ncRNA

quant_DGE_clean <-
  readRDS(file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE_edgeRfilter.RDS")

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

design <- model.matrix( ~ Passage * Treatment,
                        data = table_design)

colnames(design) <- make.names(colnames(design))

## Estimate correlation within cell pops

dupcor <- duplicateCorrelation(quant_DGE_clean$counts,
                               design,
                               block = quant_DGE_clean$samples$cell_line)

# Use ComBat-seq to correct for batch fx
# Covariates considered: Passage, Treatment

quant_DGE_batchcor_withcovariates <-
  sva::ComBat_seq(quant_DGE_clean$counts,
                  batch = quant_DGE_clean$samples$run_date,
                  covar_mod = design)

quant_DGE_clean_batchcor <- quant_DGE_clean
quant_DGE_clean_batchcor$counts <- quant_DGE_batchcor_withcovariates

# Subset data to just D3

quant_DGE_clean_batchcor_subset <- quant_DGE_clean_batchcor[, which(quant_DGE_clean_batchcor$samples$Day == "D3")]

table_design <- quant_DGE_clean_batchcor_subset$samples %>%
  mutate(Passage = factor(Passage,
                          levels = c("P5", "P7", "P13"))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Untreated", "Treated")))

design <- model.matrix( ~ Passage * Treatment,
                        data = table_design)

colnames(design) <- make.names(colnames(design))

## Apply voom transformation
### And output mean-variance trend plot

png(
  "./output/plots_QC/voom mean-variance trend.png",
  width = 20,
  height = 20,
  units = 'cm',
  res = 600
)

quant_DGE_voom <- voom(quant_DGE_clean_batchcor_subset,
                       design,
                       plot = TRUE)

dev.off()


## Duplicate correlation not necessary because cell_line is treated as a fixed (additive) effect

fit <- lmFit(quant_DGE_voom,
             design,
             block = quant_DGE_voom$targets$cell_line,
             correlation = dupcor$consensus.correlation) %>%
  eBayes()

summary(decideTests(fit))

# Define posthoc tests
## Between treatments, at each passage
## These contrasts matrix only work with the ~Passage * Treatment design with data:
## All 3 passages
## D3 only
## Hep vs Ctrl
## 20176 + 21558

contrasts_custom <- cbind(
  c(0, 1, 0, 0, 0, 0),
  c(0, 1, 0, 0, 1, 0),
  c(0, 0, 1, 0, 0, 0),
  c(0, 0, 1, 0, 0, 1),
  c(0,-1, 1, 0, 0, 0),
  c(0,-1, 1, 0,-1, 1),
  c(0, 0, 0, 1, 0, 0),
  c(0, 0, 0, 1, 1, 0),
  c(0, 0, 0, 1, 0, 1)
)

rownames(contrasts) <- colnames(design)
colnames(contrasts) <- c(
  "P7vsP5_UT_avg",
  "P7vsP5_TR_avg",
  "P13vsP7_UT_avg",
  "P13vsP7_TR_avg",
  "P13vsP5_UT_avg",
  "P13vsP5_TR_avg",
  "P5_avg",
  "P7_avg",
  "P13_avg"
)

# Test for desired contrasts
## Compare makeContrasts with the manual contrasts

fit_contrasts <- contrasts.fit(fit,
                               contrasts) %>%
  eBayes()


summary(decideTests(fit_contrasts))
