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

design <- model.matrix( ~ Passage * Treatment * cell_line,
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

quant_DGE_clean_batchcor_subset <- quant_DGE_clean_batchcor[, which(quant_DGE_clean_batchcor$samples$Day == "D3")]

table_design <- quant_DGE_clean_batchcor_subset$samples %>%
  mutate(Passage = factor(Passage,
                          levels = c("P5", "P7", "P13"))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Untreated", "Treated")))

design <- model.matrix( ~ Passage * Treatment * cell_line,
                        data = table_design)

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
## These contrasts matrix only work with the ~Passage*Treatment*cell_line design with data:
## All 3 passages
## D3 only
## Hep vs Ctrl
## 20176 + 21558

matrix_contrasts_passage <- cbind(
  c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  c(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),
  c(0, 1, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0),
  c(0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0),
  c(0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0),
  c(0, 1, 0, 0, 0, 1, 0, 0.5, 0, 0, 0.5, 0),
  c(0,-1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  c(0,-1, 1, 0, 0, 0, 0,-1, 1, 0, 0, 0),
  c(0,-1, 1, 0, 0, 0, 0,-0.5, 0.5, 0, 0, 0),
  c(0,-1, 1, 0, 0,-1, 1, 0, 0, 0, 0, 0),
  c(0,-1, 1, 0, 0,-1, 1,-1, 1, 0,-1, 1),
  c(0,-1, 1, 0, 0,-1, 1,-0.5, 0.5, 0,-0.5, 0.5),
  c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  c(0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0),
  c(0, 0, 1, 0, 0, 0, 0, 0, 0.5, 0, 0, 0),
  c(0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0),
  c(0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1),
  c(0, 0, 1, 0, 0, 0, 1, 0, 0.5, 0, 0, 0.5)
)

rownames(matrix_contrasts_passage) <- colnames(design)
colnames(matrix_contrasts_passage) <- c(
  "P7vsP5_UT_20176",
  "P7vsP5_UT_21558",
  "P7vsP5_UT_avg",
  "P7vsP5_TR_20176",
  "P7vsP5_TR_21558",
  "P7vsP5_TR_avg",
  "P13vsP7_UT_20176",
  "P13vsP7_UT_21558",
  "P13vsP7_UT_avg",
  "P13vsP7_TR_20176",
  "P13vsP7_TR_21558",
  "P13vsP7_TR_avg",
  "P13vsP5_UT_20176",
  "P13vsP5_UT_21558",
  "P13vsP5_UT_avg",
  "P13vsP5_TR_20176",
  "P13vsP5_TR_21558",
  "P13vsP5_TR_avg"
)

matrix_contrasts_treatment <- cbind(
  c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  c(0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0),
  c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0.5, 0, 0),
  c(0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0),
  c(0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0),
  c(0, 0, 0, 1, 0, 1, 0, 0, 0, 0.5, 0.5, 0),
  c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0),
  c(0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1),
  c(0, 0, 0, 1, 0, 0, 1, 0, 0, 0.5, 0, 0.5)
)

rownames(matrix_contrasts_treatment) <- colnames(design)
colnames(matrix_contrasts_treatment) <- c(
  "P5_20176",
  "P5_21558",
  "P5_avg",
  "P7_20176",
  "P7_21558",
  "P7_avg",
  "P13_20176",
  "P13_21558",
  "P13_avg"
)

fit_treatment <- contrasts.fit(fit,
                               matrix_contrasts_treatment) %>%
  eBayes()

# Test for desired contrasts
## Effect of treatment at each passage

fit_treatment <- contrasts.fit(fit,
                               matrix_contrasts_treatment) %>%
  eBayes()

## Change of expression over passages (at each treatment)

fit_passage <- contrasts.fit(fit,
                             matrix_contrasts_passage) %>%
  eBayes()

