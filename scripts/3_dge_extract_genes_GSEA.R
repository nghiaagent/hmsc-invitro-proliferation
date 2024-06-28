# Extract ordered gene lists for comparative GSEA analysis 

# These lists show expression changes between Days 3 and 5

# D5vsD3 @ Phase A UT: Coef  7
# D5vsD3 @ Phase B UT: Coef  8
# D5vsD3 @ Phase C UT: Coef  9
# D5vsD3 @ Phase A T:  Coef  10
# D5vsD3 @ Phase B T:  Coef  11
# D5vsD3 @ Phase C T:  Coef  12

# These lists show expression changes between passages at D3 and D5

#  P7vsP5 @ D3 UT: Coef  13
# P13vsP7 @ D3 UT: Coef  14
# P13vsP5 @ D3 UT: Coef  15
#  P7vsP5 @ D3  T: Coef  16
# P13vsP7 @ D3  T: Coef  17
# P13vsP5 @ D3  T: Coef  18
#  P7vsP5 @ D5 UT: Coef  19
# P13vsP7 @ D5 UT: Coef  20
# P13vsP5 @ D5 UT: Coef  21
#  P7vsP5 @ D5  T: Coef  22
# P13vsP7 @ D5  T: Coef  23
# P13vsP5 @ D5  T: Coef  24

# These lists show effect of heparin changing between passages at each day

# P7vsP5  @ D3 UT vs. T: Coefs 25
# P13vsP7 @ D3 UT vs. T: Coefs 26
# P13vsP5 @ D3 UT vs. T: Coefs 27
# P7vsP5  @ D5 UT vs. T: Coefs 28
# P13vsP7 @ D5 UT vs. T: Coefs 29
# P13vsP5 @ D5 UT vs. T: Coefs 30

# These lists show effect of heparin changing between Days 3 and 5

# D5vsD3 @ Phase A UT vs. T: Coef  37
# D5vsD3 @ Phase B UT vs. T: Coef  38
# D5vsD3 @ Phase C UT vs. T: Coef  39

# Define function to create list of lists, one with ENSEMBL ID, one with ENTREZ ID

extract_geneList <- function(fit, coef) {
  top <- topTable(fit, coef = coef, number = Inf) %>%
    arrange(desc(logFC)) %>%
    dplyr::select(c(1:8))
  
  list_GENEID <- top$logFC
  names(list_GENEID) <- as.character(top$GENEID)
  
  list_ENTREZID <- na.omit(top)$logFC
  names(list_ENTREZID) <- as.character(na.omit(top)$ENTREZID)
  
  list("GENEID" = list_GENEID, "ENTREZID" = list_ENTREZID)
}

# Ordered gene list, D5vsD3

genes_D5vsD3_UT_P5 <- extract_geneList(fit_contrasts, 7)

genes_D5vsD3_UT_P7 <- extract_geneList(fit_contrasts, 8)

genes_D5vsD3_UT_P13 <- extract_geneList(fit_contrasts, 9)

genes_D5vsD3_T_P5 <- extract_geneList(fit_contrasts, 10)

genes_D5vsD3_T_P7 <- extract_geneList(fit_contrasts, 11)

genes_D5vsD3_T_P13 <- extract_geneList(fit_contrasts, 12)

# Ordered gene list, between passages

genes_P7vsP5_UT_D3 <- extract_geneList(fit_contrasts, 13)

genes_P13vsP7_UT_D3 <- extract_geneList(fit_contrasts, 14)

genes_P13vsP5_UT_D3 <- extract_geneList(fit_contrasts, 15)

genes_P7vsP5_T_D3 <- extract_geneList(fit_contrasts, 16)

genes_P13vsP7_T_D3 <- extract_geneList(fit_contrasts, 17)

genes_P13vsP5_T_D3 <- extract_geneList(fit_contrasts, 18)

genes_P7vsP5_UT_D5 <- extract_geneList(fit_contrasts, 19)

genes_P13vsP7_UT_D5 <- extract_geneList(fit_contrasts, 20)

genes_P13vsP5_UT_D5 <- extract_geneList(fit_contrasts, 21)

genes_P7vsP5_T_D5 <- extract_geneList(fit_contrasts, 22)

genes_P13vsP7_T_D5 <- extract_geneList(fit_contrasts, 23)

genes_P13vsP5_T_D5 <- extract_geneList(fit_contrasts, 24)

# Ordered gene list, between passages at day 3

genes_P7vsP5_D3_TvsUT <- extract_geneList(fit_contrasts, 25)

genes_P13vsP7_D3_TvsUT <- extract_geneList(fit_contrasts, 26)

genes_P13vsP5_D3_TvsUT <- extract_geneList(fit_contrasts, 27)

# Ordered gene list, between passages at day 5

genes_P7vsP5_D5_TvsUT <- extract_geneList(fit_contrasts, 28)

genes_P13vsP7_D5_TvsUT <- extract_geneList(fit_contrasts, 29)

genes_P13vsP5_D5_TvsUT <- extract_geneList(fit_contrasts, 30)

# Ordered gene list, D5vsD3

genes_D5vsD3_P5_TvsUT <- extract_geneList(fit_contrasts, 37)

genes_D5vsD3_P7_TvsUT <- extract_geneList(fit_contrasts, 38)

genes_D5vsD3_P13_TvsUT <- extract_geneList(fit_contrasts, 39)
