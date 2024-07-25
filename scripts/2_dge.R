### Don't source this file by itself; call in from another file after running env_prep.R

# Load dataset

quant_DGE_clean <-
  readRDS(file = "./output/data_expression/pre_DGE/quant_cDNA_DGE_filter.RDS")

# Define design matrix for batch correction

## Treat time points and treatments as fixed effects
## Treat cell line as an additive factor

table_design <- quant_DGE_clean$samples %>%
  mutate(Passage = factor(Passage,
                          levels = c("P5", "P7", "P13"))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Untreated", "Treated")))

design <- model.matrix(~ condition_ID + cell_line,
                       data = table_design)

# Use ComBat-seq to correct for batch fx
# Batch fx corrected while considering: condition, cell population

quant_DGE_batchcor_E <-
  sva::ComBat_seq(quant_DGE_clean$counts,
                  batch = table_design$run_date,
                  covar_mod = design)

quant_DGE_clean_batchcor <- quant_DGE_clean
quant_DGE_clean_batchcor$counts <- quant_DGE_batchcor_E
rm(quant_DGE_batchcor_E)

# Define factors for design matrix

table_design <- quant_DGE_clean_batchcor$samples %>%
  mutate(Day = factor(Day,
                      levels = c("D3", "D5"))) %>%
  mutate(Passage = factor(Passage,
                          levels = c("P5", "P7", "P13"))) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("Untreated", "Treated"))) %>%
  mutate(condition_ID = factor(
    condition_ID,
    levels = c(
      "P5D3Untreated",
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
      "P13D5Treated"
    )
  ))

# Define design matrix for limma

## Treat time points and treatments (condition) as fixed effects
## Treat cell line as an additive factor
## Include batch as an additive factor

design <- model.matrix(~ condition_ID + run_date + cell_line,
                       data = table_design)

colnames(design) <- make.names(colnames(design))

# Apply voom transformation, data are now logCPM
# Output mean-variance trend plot

png(
  "./output/plots_QC/voom mean-variance trend.png",
  width = 18,
  height = 12,
  units = 'cm',
  res = 400
)

quant_DGE_voom <-
  voom(quant_DGE_clean_batchcor,
       design,
       plot = TRUE)

dev.off()

# Apply limma model fit

fit <- lmFit(quant_DGE_voom,
             design) %>%
  eBayes()

summary(decideTests(fit))

# Define contrasts
## Coefs 1 - 6: Treatment at each timepoint
## Coefs 7 - 12: Day at each timepoint x treatment
## Coefs 13 - 24: Passage at each day x treatment
## Coefs 25 - 30: Passage x Treatment at each day
## Coefs 31 - 36: Passage x Day at treatment
## Coefs 37 - 39: Day x Treatment at Passage

matrix_contrasts <- makeContrasts(
  
  ## Coefs 1 - 6: Treatment at each timepoint
  Trt_P5_D3 = condition_IDP5D3Treated - 0,
  Trt_P5_D5 = condition_IDP5D5Treated - condition_IDP5D5Untreated,
  Trt_P7_D3 = condition_IDP7D3Treated - condition_IDP7D3Untreated,
  Trt_P7_D5 = condition_IDP7D5Treated - condition_IDP7D5Untreated,
  Trt_P13_D3 = condition_IDP13D3Treated - condition_IDP13D3Untreated,
  Trt_P13_D5 = condition_IDP13D5Treated - condition_IDP13D5Untreated,
  
  ## Coefs 7 - 12: Day at each timepoint x treatment
  D5vsD3_UT_P5 = condition_IDP5D5Untreated - 0,
  D5vsD3_UT_P7 = condition_IDP7D5Untreated - condition_IDP7D3Untreated,
  D5vsD3_UT_P13 = condition_IDP13D5Untreated - condition_IDP13D3Untreated,
  D5vsD3_T_P5 = condition_IDP5D5Treated - condition_IDP5D3Treated,
  D5vsD3_T_P7 = condition_IDP7D5Treated - condition_IDP7D3Treated,
  D5vsD3_T_P13 = condition_IDP13D5Treated - condition_IDP13D3Treated,
  
  ## Coefs 13 - 24: Passage at each day x treatment
  P7vsP5_UT_D3 = condition_IDP7D3Untreated - 0,
  P13vsP7_UT_D3 = condition_IDP13D3Untreated - condition_IDP7D3Untreated,
  P13vsP5_UT_D3 = condition_IDP13D3Untreated - 0,
  P7vsP5_T_D3 = condition_IDP7D3Treated - condition_IDP5D3Treated,
  P13vsP7_T_D3 = condition_IDP13D3Treated - condition_IDP7D3Treated,
  P13vsP5_T_D3 = condition_IDP13D3Treated - condition_IDP5D3Treated,
  P7vsP5_UT_D5 = condition_IDP7D5Untreated - condition_IDP5D5Untreated,
  P13vsP7_UT_D5 = condition_IDP13D5Untreated - condition_IDP7D5Untreated,
  P13vsP5_UT_D5 = condition_IDP13D5Untreated - condition_IDP5D5Untreated,
  P7vsP5_T_D5 = condition_IDP7D5Treated - condition_IDP5D5Treated,
  P13vsP7_T_D5 = condition_IDP13D5Treated - condition_IDP7D5Treated,
  P13vsP5_T_D5 = condition_IDP13D5Treated - condition_IDP5D5Treated,
  
  ## Coefs 25 - 30: Passage x Treatment at each day
  P7vsP5_TvsUT_D3 = (condition_IDP7D3Treated - condition_IDP7D3Untreated) - (condition_IDP5D3Treated - 0),
  P13vsP7_TvsUT_D3 = (condition_IDP13D3Treated - condition_IDP13D3Untreated) - (condition_IDP7D3Treated - condition_IDP7D3Untreated),
  P13vsP5_TvsUT_D3 = (condition_IDP13D3Treated - condition_IDP13D3Untreated) - (condition_IDP5D3Treated - 0),
  P7vsP5_TvsUT_D5 = (condition_IDP7D5Treated - condition_IDP7D5Untreated) - (condition_IDP5D5Treated - condition_IDP5D5Untreated),
  P13vsP7_TvsUT_D5 = (condition_IDP13D5Treated - condition_IDP13D5Untreated) - (condition_IDP7D5Treated - condition_IDP7D5Untreated),
  P13vsP5_TvsUT_D5 = (condition_IDP13D5Treated - condition_IDP13D5Untreated) - (condition_IDP5D5Treated - condition_IDP5D5Untreated),
  
  ## Coefs 31 - 36: Passage x Day at treatment
  P7vsP5_D5vsD3_UT = (condition_IDP7D5Untreated - condition_IDP5D5Untreated) - (condition_IDP7D3Untreated - 0),
  P13vsP7_D5vsD3_UT = (condition_IDP13D5Untreated - condition_IDP7D5Untreated) - (condition_IDP13D3Untreated - condition_IDP7D3Untreated),
  P13vsP5_D5vsD3_UT = (condition_IDP13D5Untreated - condition_IDP5D5Untreated) - (condition_IDP13D3Untreated - 0),
  P7vsP5_D5vsD3_T = (condition_IDP7D5Treated - condition_IDP5D5Treated) - (condition_IDP7D3Treated - condition_IDP5D3Treated),
  P13vsP7_D5vsD3_T = (condition_IDP13D5Treated - condition_IDP7D5Treated) - (condition_IDP13D3Treated - condition_IDP7D3Treated),
  P13vsP5_D5vsD3_T = (condition_IDP13D5Treated - condition_IDP5D5Treated) - (condition_IDP13D3Treated - condition_IDP5D3Treated),
  
  ## Coefs 37 - 39: Day x Treatment at Passage
  D5vsD3_TvsUT_P5 = (condition_IDP5D5Treated - condition_IDP5D5Untreated) - (condition_IDP5D3Treated - 0),
  D5vsD3_TvsUT_P7 = (condition_IDP7D5Treated - condition_IDP7D5Untreated) - (condition_IDP7D3Treated - condition_IDP7D3Untreated),
  D5vsD3_TvsUT_P13 = (condition_IDP13D5Treated - condition_IDP13D5Untreated) - (condition_IDP13D3Treated - condition_IDP13D3Untreated),
  
  levels = design
)

# Test for desired contrasts

fit_contrasts <- contrasts.fit(fit,
                               matrix_contrasts) %>%
  eBayes()

# Save data

saveRDS(quant_DGE_voom,
        file = "./output/data_expression/post_DGE/quant_DGE_voom.RDS")

saveRDS(fit_contrasts,
        file = "./output/data_expression/post_DGE/fit_contrasts.RDS")

saveRDS(fit,
        file = "./output/data_expression/post_DGE/fit.RDS")