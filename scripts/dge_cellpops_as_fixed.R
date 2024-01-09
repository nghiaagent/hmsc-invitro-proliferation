### Don't source this file by itself; call in from another file after running env_prep.R
# Load dataset
## Dataset has both mRNA and ncRNA

quant_DGE_clean <-
  readRDS(file = "./output/quant_cDNA_ncRNA_ENSEMBL_DGE_filter.RDS")

## MDS / PCA plot

Glimma::glMDSPlot(quant_DGE_clean,
                  labels = quant_DGE_clean$samples$Sample_ID,
                  groups = quant_DGE_clean$samples)

### Make design matrix

## Treat cell line as an additive factor
## Treat time points and treatments as fixed effects

table_design <- quant_DGE_clean$samples %>%
  mutate(Passage = factor(Passage,
                          levels = c("P5", "P7", "P13"))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Untreated", "Treated")))

design <- model.matrix( ~ condition_ID,
                        data = table_design)

# Use ComBat-seq to correct for batch fx
# Covariates considered: condition, cell population

quant_DGE_batchcor_withcovariates <-
  sva::ComBat_seq(quant_DGE_clean$counts,
                  batch = table_design$run_date,
                  covar_mod = design)

quant_DGE_clean_batchcor <- quant_DGE_clean
quant_DGE_clean_batchcor$counts <- quant_DGE_batchcor_withcovariates

quant_DGE_clean_batchcor_subset <- quant_DGE_clean_batchcor
  
table_design <- quant_DGE_clean_batchcor_subset$samples %>%
  mutate(Day = factor(Day,
                      levels = c("D3", "D5"))) %>%
  mutate(Passage = factor(Passage,
                          levels = c("P5", "P7", "P13"))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Untreated", "Treated"))) %>%
  mutate(condition_ID = factor(condition_ID,
                               levels = c("P5D3Untreated",
                                          "P5D3Treated",
                                          "P5D5Untreated",
                                          "P5D5Treated",
                                          "P7D3Untreated",
                                          "P7D3Treated",
                                          "P7D5Untreated",
                                          "P7D5Treated",
                                          "P13D3Untreated",
                                          "P13D3Treated",
                                          "P13D5Untreated",
                                          "P13D5Treated")))

design <- model.matrix( ~ condition_ID + run_date,
                        data = table_design)

## Apply voom transformation
### And output mean-variance trend plot

png(
  "./output/plots_QC/voom mean-variance trend.png",
  width = 60,
  height = 20,
  units = 'cm',
  res = 400
)


quant_DGE_voom <-
  voom(quant_DGE_clean_batchcor_subset,
                         design,
                         plot = TRUE)

dev.off()

## Apply limma model fit

fit <- lmFit(quant_DGE_voom,
             design) %>%
  eBayes()

summary(decideTests(fit))

# Define posthoc tests
## Between treatments, at each passage
## Each coefficient calculates the difference between that condition and the baseline (P5D3Untreated)

matrix_contrasts <- makeContrasts(
  Trt_P5_D3 = condition_IDP5D3Treated - 0,
  Trt_P7_D3 = condition_IDP7D3Treated - condition_IDP7D3Untreated,
  Trt_P13_D3 = condition_IDP13D3Treated - condition_IDP13D3Untreated,
  P7vsP5_UT_D3 = condition_IDP7D3Untreated - 0,
  P13vsP5_UT_D3 = condition_IDP13D3Untreated - 0,
  P13vsP7_UT_D3 = condition_IDP13D3Untreated - condition_IDP7D3Untreated,
  P7vsP5_T_D3 = condition_IDP7D3Treated - condition_IDP5D3Treated,
  P13vsP5_T_D3 = condition_IDP13D3Treated - condition_IDP5D3Treated,
  P13vsP7_T_D3 = condition_IDP13D3Treated - condition_IDP7D3Treated,
  Trt_P5_D5 = condition_IDP5D5Treated - condition_IDP5D5Untreated,
  Trt_P7_D5 = condition_IDP7D5Treated - condition_IDP7D5Untreated,
  Trt_P13_D5 = condition_IDP13D5Treated - condition_IDP13D5Untreated,
  P7vsP5_UT_D5 = condition_IDP7D5Untreated - condition_IDP5D5Untreated,
  P13vsP5_UT_D5 = condition_IDP13D5Untreated - condition_IDP5D5Untreated,
  P13vsP7_UT_D5 = condition_IDP13D5Untreated - condition_IDP7D5Untreated,
  P7vsP5_T_D5 = condition_IDP7D5Treated - condition_IDP5D5Treated,
  P13vsP5_T_D5 = condition_IDP13D5Treated - condition_IDP5D5Treated,
  P13vsP7_T_D5 = condition_IDP13D5Treated - condition_IDP7D5Treated,
  levels = design
)


# Test for desired contrasts

fit_contrasts <- contrasts.fit(fit,
                               matrix_contrasts) %>%
  eBayes()

summary(decideTests(fit_contrasts))

