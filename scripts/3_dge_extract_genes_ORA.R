# Extract gene names from limma fit for ORA
# Between passages at D3
# Venn diagrams
# D5vsD3 @ Phase A UT vs. D5vsD3 @ Phase A T: Coefs  7; 10
# D5vsD3 @ Phase B UT vs. D5vsD3 @ Phase B T: Coefs  8; 11
# D5vsD3 @ Phase C UT vs. D5vsD3 @ Phase C T: Coefs  9; 12
# P7vsP5 @ D3 UT vs. P7vsP5 @ D3 T:           Coefs 13; 16
# P13vsP7 @ D3 UT vs. P7vsP5 @ D3 T:          Coefs 14; 17
# P13vsP5 @ D3 UT vs. P7vsP5 @ D3 T:          Coefs 15; 18
# P7vsP5 @ D5 UT vs. P7vsP5 @ D5 T:           Coefs 19; 22
# P13vsP7 @ D5 UT vs. P7vsP5 @ D5 T:          Coefs 20; 23
# P13vsP5 @ D5 UT vs. P7vsP5 @ D5 T:          Coefs 21; 24

tests <- decideTests(fit_contrasts)

# Between passages, untreated cells

genes_P7vsP5_D3_Untreated <- 
  quant_DGE_voom$genes[which(tests[, 13] != 0),]

genes_P13vsP7_D3_Untreated <- 
  quant_DGE_voom$genes[which(tests[, 14] != 0),]

genes_P13vsP5_D3_Untreated <- 
  quant_DGE_voom$genes[which(tests[, 15] != 0),]

# These Venn diagrams show overlap in changes between passages in untreated cells

## D3

genes_D3_Untreated_P7vsP5_only <-
  quant_DGE_voom$genes[which((tests[, 13] == 1 &
                                tests[, 14] != 1 &
                                tests[, 15] != 1) |
                               (tests[, 13] == -1 &
                                  tests[, 14] != -1 &
                                  tests[, 15] != -1)),]

genes_D3_Untreated_P13vsP7_only <-
  quant_DGE_voom$genes[which((tests[, 13] != 1 &
                                tests[, 14] == 1 &
                                tests[, 15] != 1) |
                               (tests[, 13] != -1 &
                                  tests[, 14] == -1 &
                                  tests[, 15] != -1)),]

genes_D3_Untreated_P13vsP5_only <-
  quant_DGE_voom$genes[which((tests[, 13] != 1 &
                                tests[, 14] != 1 &
                                tests[, 15] == 1) |
                               (tests[, 13] != -1 &
                                  tests[, 14] != -1 &
                                  tests[, 15] == -1)),]

genes_D3_Untreated_P13vsP5_and_P13vsP7_only <-
  quant_DGE_voom$genes[which((tests[, 13] != 1 &
                                tests[, 14] == 1 &
                                tests[, 15] == 1) |
                               (tests[, 13] != -1 &
                                  tests[, 14] == -1 &
                                  tests[, 15] == -1)),]


## D5

genes_D5_Untreated_P7vsP5_only <-
  quant_DGE_voom$genes[which((tests[, 19] == 1 &
                                tests[, 20] != 1 &
                                tests[, 21] != 1) |
                               (tests[, 19] == -1 &
                                  tests[, 20] != -1 &
                                  tests[, 21] != -1)),]

genes_D5_Untreated_P13vsP7_only <-
  quant_DGE_voom$genes[which((tests[, 19] != 1 &
                                tests[, 20] == 1 &
                                tests[, 21] != 1) |
                               (tests[, 19] != -1 &
                                  tests[, 20] == -1 &
                                  tests[, 21] != -1)),]

genes_D5_Untreated_P13vsP5_only <-
  quant_DGE_voom$genes[which((tests[, 19] != 1 &
                                tests[, 20] != 1 &
                                tests[, 21] == 1) |
                               (tests[, 19] != -1 &
                                  tests[, 20] != -1 &
                                  tests[, 21] == -1)),]

genes_D5_Untreated_P13vsP5_and_P13vsP7_only <-
  quant_DGE_voom$genes[which((tests[, 19] != 1 &
                                tests[, 20] == 1 &
                                tests[, 21] == 1) |
                               (tests[, 19] != -1 &
                                  tests[, 20] == -1 &
                                  tests[, 21] == -1)),]


# Export data

## Difference between passages

### P5vsP7

fwrite(
  genes_P7vsP5_D3_Untreated$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_P7vsP5_D3_Untreated_ENSEMBL.txt",
  na = ''
)

fwrite(
  genes_P7vsP5_D3_Untreated$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_P7vsP5_D3_Untreated_GENENAME.txt",
  na = ''
)

### P7vsP13

fwrite(
  genes_P13vsP7_D3_Untreated$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_P13vsP7_D3_Untreated_ENSEMBL.txt",
  na = ''
)

fwrite(
  genes_P13vsP7_D3_Untreated$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_P13vsP7_D3_Untreated_GENENAME.txt",
  na = ''
)

### P5vsP13

fwrite(
  genes_P13vsP5_D3_Untreated$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_P13vsP5_D3_Untreated_ENSEMBL.txt",
  na = ''
)

fwrite(
  genes_P13vsP5_D3_Untreated$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_P13vsP5_D3_Untreated_GENENAME.txt",
  na = ''
)

## Overlap between passages in untreated cells, D3

fwrite(
  genes_D3_Untreated_P13vsP5_and_P13vsP7_only$GENEID %>% list(),
  file = "./output/list_genes/venn_P13vsP7vsP5_UT_D3/P13vsP5_and_P13vsP7_only_ENSEMBL.txt",
  na = ''
)

fwrite(
  genes_D3_Untreated_P13vsP5_and_P13vsP7_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP7vsP5_UT_D3/P13vsP5_and_P13vsP7_only_GENENAME.txt",
  na = ''
)

fwrite(genes_D3_Untreated_P13vsP5_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7vsP5_UT_D3/P13vsP5_only_ENSEMBL.txt",
       na = '')

fwrite(
  genes_D3_Untreated_P13vsP5_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP7vsP5_UT_D3/P13vsP5_only_GENENAME.txt",
  na = ''
)

fwrite(genes_D3_Untreated_P13vsP7_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7vsP5_UT_D3/P13vsP7_only_ENSEMBL.txt",
       na = '')

fwrite(
  genes_D3_Untreated_P13vsP7_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP7vsP5_UT_D3/P13vsP7_only_GENENAME.txt",
  na = ''
)

fwrite(genes_D3_Untreated_P7vsP5_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7vsP5_UT_D3/P7vsP5_only_ENSEMBL.txt",
       na = '')

fwrite(genes_D3_Untreated_P7vsP5_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP7vsP5_UT_D3/P7vsP5_only_GENENAME.txt",
       na = '')

## Overlap between passages in untreated cells, D5

fwrite(
  genes_D5_Untreated_P13vsP5_and_P13vsP7_only$GENEID %>% list(),
  file = "./output/list_genes/venn_P13vsP7vsP5_UT_D5/P13vsP5_and_P13vsP7_only_ENSEMBL.txt",
  na = ''
)

fwrite(
  genes_D5_Untreated_P13vsP5_and_P13vsP7_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP7vsP5_UT_D5/P13vsP5_and_P13vsP7_only_GENENAME.txt",
  na = ''
)

fwrite(genes_D5_Untreated_P13vsP5_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7vsP5_UT_D5/P13vsP5_only_ENSEMBL.txt",
       na = '')

fwrite(
  genes_D5_Untreated_P13vsP5_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP7vsP5_UT_D5/P13vsP5_only_GENENAME.txt",
  na = ''
)

fwrite(genes_D5_Untreated_P13vsP7_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7vsP5_UT_D5/P13vsP7_only_ENSEMBL.txt",
       na = '')

fwrite(
  genes_D5_Untreated_P13vsP7_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP7vsP5_UT_D5/P13vsP7_only_GENENAME.txt",
  na = ''
)

fwrite(genes_D5_Untreated_P7vsP5_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7vsP5_UT_D5/P7vsP5_only_ENSEMBL.txt",
       na = '')

fwrite(genes_D5_Untreated_P7vsP5_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP7vsP5_UT_D5/P7vsP5_only_GENENAME.txt",
       na = '')
