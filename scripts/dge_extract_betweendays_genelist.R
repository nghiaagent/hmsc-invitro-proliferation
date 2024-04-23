# Extract ordered gene list for comparative GSEA analysis 

# These lists show expression changes between Days 3 and 5

# D5vsD3 @ Phase A UT: Coef  7
# D5vsD3 @ Phase B UT: Coef  8
# D5vsD3 @ Phase C UT: Coef  9
# D5vsD3 @ Phase A T:  Coef  10
# D5vsD3 @ Phase B T:  Coef  11
# D5vsD3 @ Phase C T:  Coef  12

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

genes_D5vsD3_UT_P5 <- extract_geneList(fit_contrasts, 7)

genes_D5vsD3_UT_P7 <- extract_geneList(fit_contrasts, 8)

genes_D5vsD3_UT_P13 <- extract_geneList(fit_contrasts, 9)

genes_D5vsD3_T_P5 <- extract_geneList(fit_contrasts, 10)

genes_D5vsD3_T_P7 <- extract_geneList(fit_contrasts, 11)

genes_D5vsD3_T_P13 <- extract_geneList(fit_contrasts, 12)

