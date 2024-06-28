# Extract ordered gene lists for camera 
# These lists show expression changes between passages at D3 and D5

#  P7vsP5 @ D3 UT: Coef  13
# P13vsP7 @ D3 UT: Coef  14
# P13vsP5 @ D3 UT: Coef  15

# These lists show effect of heparin changing between passages at day 3

# P7vsP5  @ D3 UT vs. T: Coefs 25
# P13vsP7 @ D3 UT vs. T: Coefs 26
# P13vsP5 @ D3 UT vs. T: Coefs 27

# These lists show effect of heparin at each passage

# P5 D3: Coef  1
# P7 D3: Coef  3
# P13D3: Coef  5

extract_genelist <- function(fit, coef) {
  top <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  
  list_GENEID <- top$t
  names(list_GENEID) <- as.character(top$GENEID)
  return(list_GENEID)
}

# Ordered gene list, between passages

genes_P7vsP5_UT_D3 <- extract_genelist(fit_contrasts, 13)

genes_P13vsP7_UT_D3 <- extract_genelist(fit_contrasts, 14)

genes_P13vsP5_UT_D3 <- extract_genelist(fit_contrasts, 15)


# Ordered gene list, between passages at day 3

genes_P7vsP5_D3_TvsUT <- extract_genelist(fit_contrasts, 25)

genes_P13vsP7_D3_TvsUT <- extract_genelist(fit_contrasts, 26)

genes_P13vsP5_D3_TvsUT <- extract_genelist(fit_contrasts, 27)

# Ranked gene list, between treatments at day 3

genes_TvsUT_P5_D3 <- extract_genelist(fit_contrasts, 1)

genes_TvsUT_P7_D3 <- extract_genelist(fit_contrasts, 3)

genes_TvsUT_P13_D3 <- extract_genelist(fit_contrasts, 5)