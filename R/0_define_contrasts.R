here::i_am("R/0_define_contrasts.R")

# Define contrasts for DESeq2
list_contrasts_deseq2 <- list(
  ## Coefs 1 - 6: Treatment at each timepoint
  Trt_P5_D3 = c("condition_ID", "P5D3Treated", "P5D3Untreated"),
  Trt_P5_D5 = c("condition_ID", "P5D5Treated", "P5D5Untreated"),
  Trt_P7_D3 = c("condition_ID", "P7D3Treated", "P7D3Untreated"),
  Trt_P7_D5 = c("condition_ID", "P7D5Treated", "P7D5Untreated"),
  Trt_P13_D3 = c("condition_ID", "P13D3Treated", "P13D3Untreated"),
  Trt_P13_D5 = c("condition_ID", "P13D5Treated", "P13D5Untreated"),

  ## Coefs 7 - 12: Day at each timepoint x treatment
  D5vsD3_UT_P5 = c("condition_ID", "P5D5Untreated", "P5D3Untreated"),
  D5vsD3_UT_P7 = c("condition_ID", "P7D5Untreated", "P7D3Untreated"),
  D5vsD3_UT_P13 = c("condition_ID", "P13D5Untreated", "P13D3Untreated"),
  D5vsD3_T_P5 = c("condition_ID", "P5D5Treated", "P5D3Treated"),
  D5vsD3_T_P7 = c("condition_ID", "P7D5Treated", "P7D3Treated"),
  D5vsD3_T_P13 = c("condition_ID", "P13D5Treated", "P13D3Treated"),

  ## Coefs 13 - 24: Passage at each day x treatment
  P7vsP5_UT_D3 = c("condition_ID", "P7D3Untreated", "P5D3Untreated"),
  P13vsP7_UT_D3 = c("condition_ID", "P13D3Untreated", "P7D3Untreated"),
  P13vsP5_UT_D3 = c("condition_ID", "P13D3Untreated", "P5D3Untreated"),
  P7vsP5_T_D3 = c("condition_ID", "P7D3Treated", "P5D3Treated"),
  P13vsP7_T_D3 = c("condition_ID", "P13D3Treated", "P7D3Treated"),
  P13vsP5_T_D3 = c("condition_ID", "P13D3Treated", "P5D3Treated"),
  P7vsP5_UT_D5 = c("condition_ID", "P7D5Untreated", "P5D5Untreated"),
  P13vsP7_UT_D5 = c("condition_ID", "P13D5Untreated", "P7D5Untreated"),
  P13vsP5_UT_D5 = c("condition_ID", "P13D5Untreated", "P5D5Untreated"),
  P7vsP5_T_D5 = c("condition_ID", "P7D5Treated", "P5D5Treated"),
  P13vsP7_T_D5 = c("condition_ID", "P13D5Treated", "P7D5Treated"),
  P13vsP5_T_D5 = c("condition_ID", "P13D5Treated", "P5D5Treated")
)
