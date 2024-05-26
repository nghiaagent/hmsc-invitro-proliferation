# Extract ordered gene list for comparative GSEA analysis 

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

# Create list of lists, one with ENSEMBL ID, one with ENTREZ ID

extract_geneList <- function(fit, coef) {
  
  top <- topTable(fit,
                   coef = coef,
                  number = Inf) %>%
    arrange(desc(logFC)) %>%
    dplyr::select(c(1:7))
  
  list_GENEID <- top$logFC
  names(list_GENEID) <- as.character(top$GENEID)
  
  list_ENTREZID <- na.omit(top)$logFC
  names(list_ENTREZID) <- as.character(na.omit(top)$ENTREZID)

  list(
    "GENEID" = list_GENEID,
    "ENTREZID" = list_ENTREZID
  )
}

# Ordered gene list, D5vsD3

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
