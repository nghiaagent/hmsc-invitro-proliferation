#### Load dataset ####

## This dataset is already filtered and normalised.

quant_DGE_clean <-
  readRDS(file = "./output/data_expression/pre_DGE/quant_cDNA_DGE_filter.RDS")

#### Define design matrix ####
## The design matrix defines our experimental design. There are multiple "valid" options.
## The selected design facilitates easy pairwise comparisons between conditions of interest.
## See [TODO: Insert model matrix guide] for alternative designs.

design <- model.matrix(~ condition_ID,
                       data = quant_DGE_clean$samples)

colnames(design) %<>% make.names(.)

#### Define contrast matrix ####
## The contrast matrix defines our desired comparisons.
## NB: 0 represents the "baseline" condition (WT, untreated).

matrix_contrasts <- makeContrasts(
  
  # Contrast  1: Change in baseline expression between KO and WT cells
  KO_vs_WT_UT          = condition_IDATP7B.KO_Untreated   - 0,
  
  # Contrasts 2 - 6: Effect of treatments in wild-type WT cells
  CU_vs_UT_WT          =  condition_IDWT_Cu                - 0,
  DPEN_vs_CU_WT        =  condition_IDWT_Cu_D.penicilamine - condition_IDWT_Cu,
  TRIEN_vs_Cu_WT       =  condition_IDWT_Cu_trientine      - condition_IDWT_Cu,
  DPEN_vs_UT_WT        =  condition_IDWT_Cu_D.penicilamine - 0,
  TRIEN_vs_UT_WT       =  condition_IDWT_Cu_trientine      - 0,
  
  # Contrasts 7  - 11: Effect of treatments in ATP7B-KO KO cells
  CU_vs_UT_KO          =  condition_IDATP7B.KO_Cu                - condition_IDATP7B.KO_Untreated,
  DPEN_vs_CU_KO        =  condition_IDATP7B.KO_Cu_D.penicilamine - condition_IDATP7B.KO_Cu,
  TRIEN_vs_CU_KO       =  condition_IDATP7B.KO_Cu_trientine      - condition_IDATP7B.KO_Cu,
  DPEN_vs_UT_KO        =  condition_IDATP7B.KO_Cu_D.penicilamine - condition_IDATP7B.KO_Untreated,
  TRIEN_vs_UT_KO       =  condition_IDATP7B.KO_Cu_trientine      - condition_IDATP7B.KO_Untreated,
  
  # Contrasts 12 - 16: Change in the effect of treatment between KO and WT cells
  CU_vs_UT_KO_vs_WT    = (condition_IDATP7B.KO_Cu                - condition_IDATP7B.KO_Untreated) - (condition_IDWT_Cu                - 0),
  DPEN_vs_CU_KO_vs_WT  = (condition_IDATP7B.KO_Cu_D.penicilamine - condition_IDATP7B.KO_Cu)        - (condition_IDWT_Cu_D.penicilamine - condition_IDWT_Cu),
  TRIEN_vs_CU_KO_vs_WT = (condition_IDATP7B.KO_Cu_trientine      - condition_IDATP7B.KO_Cu)        - (condition_IDWT_Cu_trientine      - condition_IDWT_Cu),
  DPEN_vs_UT_KO_vs_WT  = (condition_IDATP7B.KO_Cu_D.penicilamine - condition_IDATP7B.KO_Untreated) - (condition_IDWT_Cu_D.penicilamine - 0),
  TRIEN_vs_UT_KO_vs_WT = (condition_IDATP7B.KO_Cu_trientine      - condition_IDATP7B.KO_Untreated) - (condition_IDWT_Cu_trientine      - 0),
  
  levels = design
)

#### Apply voom and limma model fit ####
# Voom converts counts to logCPM (equiv. to log-expression) and
# estimates variance to allow application of linear models on heteroscedastic data.
# Limma applies a linear model to each gene to estimate its expression at each condition.
# Statistics are shared between the ~15000 linear models to increase statistical power.
# Extra reading: Google "STAT555 lecture".
# Output mean-variance trend plot


png(
  "./output/plots_QC/voom mean-variance trend.png",
  width = 18,
  height = 12,
  units = 'cm',
  res = 400
)

fit_contrasts <- voomLmFit(quant_DGE_clean,
                           design,
                           keep.EList = TRUE,
                           plot = TRUE) %>%
  eBayes() %>%
  contrasts.fit(matrix_contrasts) %>%
  eBayes() %>%
  ## Subset for protein-coding genes
  .[which(.$genes$GENETYPE == "protein-coding"), ]

fit_contrasts[["EList"]] %<>% .[which(.$genes$GENETYPE == "protein-coding"), ]

dev.off()

#### Save output ####

saveRDS(fit_contrasts, file = "./output/data_expression/post_DGE/fit_contrasts.RDS")