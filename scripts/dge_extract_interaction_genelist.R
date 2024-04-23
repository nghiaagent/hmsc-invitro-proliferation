# Extract ordered gene list from interaction terms names from Venn diagram

# These lists show effect of heparin changing between Days 3 and 5

# D5vsD3 @ Phase A UT vs. T: Coef  37
# D5vsD3 @ Phase B UT vs. T: Coef  38
# D5vsD3 @ Phase C UT vs. T: Coef  39

# These lists show effect of heparin changing between passages at each day

# P7vsP5  @ D3 UT vs. T: Coefs 25
# P13vsP7 @ D3 UT vs. T: Coefs 26
# P13vsP5 @ D3 UT vs. T: Coefs 27
# P7vsP5  @ D5 UT vs. T: Coefs 28
# P13vsP7 @ D5 UT vs. T: Coefs 29
# P13vsP5 @ D5 UT vs. T: Coefs 30

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

genes_D5vsD3_P5_TvsUT <- extract_geneList(fit_contrasts, 37)

genes_D5vsD3_P7_TvsUT <- extract_geneList(fit_contrasts, 38)

genes_D5vsD3_P13_TvsUT <- extract_geneList(fit_contrasts, 39)

# Ordered gene list, between passages at day 3

genes_P7vsP5_D3_TvsUT <- extract_geneList(fit_contrasts, 25)

genes_P13vsP7_D3_TvsUT <- extract_geneList(fit_contrasts, 26)

genes_P13vsP5_D3_TvsUT <- extract_geneList(fit_contrasts, 27)

# Ordered gene list, between passages at day 5

genes_P7vsP5_D5_TvsUT <- extract_geneList(fit_contrasts, 28)

genes_P13vsP7_D5_TvsUT <- extract_geneList(fit_contrasts, 29)

genes_P13vsP5_D5_TvsUT <- extract_geneList(fit_contrasts, 30)