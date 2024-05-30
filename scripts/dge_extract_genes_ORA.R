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

tests <- decideTests(fit_contrasts[!is.na(fit_contrasts$genes$ENTREZID),])

# Between passages, untreated cells

genes_P7vsP5_D3_Untreated <- 
  quant_DGE_voom$genes[which(tests[, 13] != 0),]

genes_P13vsP7_D3_Untreated <- 
  quant_DGE_voom$genes[which(tests[, 14] != 0),]

genes_P13vsP5_D3_Untreated <- 
  quant_DGE_voom$genes[which(tests[, 15] != 0),]

# D5vsD3 at passages

genes_D5vsD3_P5_Treated_only <-
  quant_DGE_voom$genes[which((tests[, 7] != 1 &
                                tests[, 10] == 1) |
                               (tests[, 7] != -1 &
                                  tests[, 10] == -1)),]

genes_D5vsD3_P5_Untreated_only <-
  quant_DGE_voom$genes[which((tests[, 7] == 1 &
                                tests[, 10] != 1) |
                               (tests[, 7] == -1 &
                                  tests[, 10] != -1)),]

genes_D5vsD3_P5_Overlap <-
  quant_DGE_voom$genes[which((tests[, 7] == 1 &
                                tests[, 10] == 1) |
                               (tests[, 7] == -1 &
                                  tests[, 10] == -1)),]

genes_D5vsD3_P7_Treated_only <-
  quant_DGE_voom$genes[which((tests[, 8] != 1 &
                                tests[, 11] == 1) |
                               (tests[, 8] != -1 &
                                  tests[, 11] == -1)),]

genes_D5vsD3_P7_Untreated_only <-
  quant_DGE_voom$genes[which((tests[, 8] == 1 &
                                tests[, 11] != 1) |
                               (tests[, 8] == -1 &
                                  tests[, 11] != -1)),]

genes_D5vsD3_P7_Overlap <-
  quant_DGE_voom$genes[which((tests[, 8] == 1 &
                                tests[, 11] == 1) |
                               (tests[, 8] == -1 &
                                  tests[, 11] == -1)),]


genes_D5vsD3_P13_Treated_only <-
  quant_DGE_voom$genes[which((tests[, 9] != 1 &
                                tests[, 12] == 1) |
                               (tests[, 9] != -1 &
                                  tests[, 12] == -1)),]

genes_D5vsD3_P13_Untreated_only <-
  quant_DGE_voom$genes[which((tests[, 9] == 1 &
                                tests[, 12] != 1) |
                               (tests[, 9] == -1 &
                                  tests[, 12] != -1)),]

genes_D5vsD3_P13_Overlap <-
  quant_DGE_voom$genes[which((tests[, 9] == 1 &
                                tests[, 12] == 1) |
                               (tests[, 9] == -1 &
                                  tests[, 12] == -1)),]

# Between passages at D3


genes_P7vsP5_D3_Treated_only <-
  quant_DGE_voom$genes[which((tests[, 13] != 1 &
                                tests[, 16] == 1) |
                               (tests[, 13] != -1 &
                                  tests[, 16] == -1)),]

genes_P7vsP5_D3_Untreated_only <-
  quant_DGE_voom$genes[which((tests[, 13] == 1 &
                                tests[, 16] != 1) |
                               (tests[, 13] == -1 &
                                  tests[, 16] != -1)),]

genes_P7vsP5_D3_Overlap <-
  quant_DGE_voom$genes[which((tests[, 13] == 1 &
                                tests[, 16] == 1) |
                               (tests[, 13] == -1 &
                                  tests[, 16] == -1)),]


genes_P13vsP7_D3_Treated_only <-
  quant_DGE_voom$genes[which((tests[, 14] != 1 &
                                tests[, 17] == 1) |
                               (tests[, 14] != -1 &
                                  tests[, 17] == -1)),]

genes_P13vsP7_D3_Untreated_only <-
  quant_DGE_voom$genes[which((tests[, 14] == 1 &
                                tests[, 17] != 1) |
                               (tests[, 14] == -1 &
                                  tests[, 17] != -1)),]


genes_P13vsP7_D3_Overlap <-
  quant_DGE_voom$genes[which((tests[, 14] == 1 &
                                tests[, 17] == 1) |
                               (tests[, 14] == -1 &
                                  tests[, 17] == -1)),]


genes_P13vsP5_D3_Treated_only <-
  quant_DGE_voom$genes[which((tests[, 15] != 1 &
                                tests[, 18] == 1) |
                               (tests[, 15] != -1 &
                                  tests[, 18] == -1)),]

genes_P13vsP5_D3_Untreated_only <-
  quant_DGE_voom$genes[which((tests[, 15] == 1 &
                                tests[, 18] != 1) |
                               (tests[, 15] == -1 &
                                  tests[, 18] != -1)),]

genes_P13vsP5_D3_Overlap <-
  quant_DGE_voom$genes[which((tests[, 15] == 1 &
                                tests[, 18] == 1) |
                               (tests[, 15] == -1 &
                                  tests[, 18] == -1)),]

# Between passages at D5

genes_P7vsP5_D5_Treated_only <-
  quant_DGE_voom$genes[which((tests[, 19] != 1 &
                                tests[, 22] == 1) |
                               (tests[, 19] != -1 &
                                  tests[, 22] == -1)),]

genes_P7vsP5_D5_Untreated_only <-
  quant_DGE_voom$genes[which((tests[, 19] == 1 &
                                tests[, 22] != 1) |
                               (tests[, 19] == -1 &
                                  tests[, 22] != -1)),]

genes_P7vsP5_D5_Overlap <-
  quant_DGE_voom$genes[which((tests[, 19] == 1 &
                                tests[, 22] == 1) |
                               (tests[, 19] == -1 &
                                  tests[, 22] == -1)),]


genes_P13vsP7_D5_Treated_only <-
  quant_DGE_voom$genes[which((tests[, 20] != 1 &
                                tests[, 23] == 1) |
                               (tests[, 20] != -1 &
                                  tests[, 23] == -1)),]

genes_P13vsP7_D5_Untreated_only <-
  quant_DGE_voom$genes[which((tests[, 20] == 1 &
                                tests[, 23] != 1) |
                               (tests[, 20] == -1 &
                                  tests[, 23] != -1)),]

genes_P13vsP7_D5_Overlap <-
  quant_DGE_voom$genes[which((tests[, 20] == 1 &
                                tests[, 23] == 1) |
                               (tests[, 20] == -1 &
                                  tests[, 23] == -1)),]


genes_P13vsP5_D5_Treated_only <-
  quant_DGE_voom$genes[which((tests[, 21] != 1 &
                                tests[, 24] == 1) |
                               (tests[, 21] != -1 &
                                  tests[, 24] == -1)),]

genes_P13vsP5_D5_Untreated_only <-
  quant_DGE_voom$genes[which((tests[, 21] == 1 &
                                tests[, 24] != 1) |
                               (tests[, 21] == -1 &
                                  tests[, 24] != -1)),]

genes_P13vsP5_D5_Overlap <-
  quant_DGE_voom$genes[which((tests[, 21] == 1 &
                                tests[, 24] == 1) |
                               (tests[, 21] == -1 &
                                  tests[, 24] == -1)),]


# These Venn diagrams show overlap between passages in untreated cells

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

## Between passages D3

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

## Between passages D5

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

## Between days at each passage, TvsUT

fwrite(genes_D5vsD3_P5_Untreated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P5/D5vsD3_P5_Untreated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_D5vsD3_P5_Untreated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P5/D5vsD3_P5_Untreated_only_GENENAME.txt",
       na = '')

fwrite(genes_D5vsD3_P5_Treated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P5/D5vsD3_P5_Treated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_D5vsD3_P5_Treated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P5/D5vsD3_P5_Treated_only_GENENAME.txt",
       na = '')

fwrite(genes_D5vsD3_P5_Overlap$GENEID %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P5/D5vsD3_P5_Overlap_ENSEMBL.txt",
       na = '')

fwrite(genes_D5vsD3_P5_Overlap$GENENAME %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P5/D5vsD3_P5_Overlap_GENENAME.txt",
       na = '')

fwrite(genes_D5vsD3_P7_Untreated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P7/D5vsD3_P7_Untreated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_D5vsD3_P7_Untreated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P7/D5vsD3_P7_Untreated_only_GENENAME.txt",
       na = '')

fwrite(genes_D5vsD3_P7_Treated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P7/D5vsD3_P7_Treated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_D5vsD3_P7_Treated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P7/D5vsD3_P7_Treated_only_GENENAME.txt",
       na = '')

fwrite(genes_D5vsD3_P7_Overlap$GENEID %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P7/D5vsD3_P7_Overlap_ENSEMBL.txt",
       na = '')

fwrite(genes_D5vsD3_P7_Overlap$GENENAME %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P7/D5vsD3_P7_Overlap_GENENAME.txt",
       na = '')

fwrite(genes_D5vsD3_P13_Untreated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P13/D5vsD3_P13_Untreated_only_ENSEMBL.txt",
       na = '')

fwrite(
  genes_D5vsD3_P13_Untreated_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_D5vsD3_P13/D5vsD3_P13_Untreated_only_GENENAME.txt",
  na = ''
)

fwrite(genes_D5vsD3_P13_Treated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P13/D5vsD3_P13_Treated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_D5vsD3_P13_Treated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P13/D5vsD3_P13_Treated_only_GENENAME.txt",
       na = '')

fwrite(genes_D5vsD3_P13_Overlap$GENEID %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P13/D5vsD3_P13_Overlap_ENSEMBL.txt",
       na = '')

fwrite(genes_D5vsD3_P13_Overlap$GENENAME %>% list(),
       file = "./output/list_genes/venn_D5vsD3_P13/D5vsD3_P13_Overlap_GENENAME.txt",
       na = '')

## Between passages at day 3, TvsUT

fwrite(genes_P7vsP5_D3_Untreated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D3/P7vsP5_D3_Untreated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_P7vsP5_D3_Untreated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D3/P7vsP5_D3_Untreated_only_GENENAME.txt",
       na = '')

fwrite(genes_P7vsP5_D3_Treated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D3/P7vsP5_D3_Treated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_P7vsP5_D3_Treated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D3/P7vsP5_D3_Treated_only_GENENAME.txt",
       na = '')

fwrite(genes_P7vsP5_D3_Overlap$GENEID %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D3/P7vsP5_D3_Overlap_ENSEMBL.txt",
       na = '')

fwrite(genes_P7vsP5_D3_Overlap$GENENAME %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D3/P7vsP5_D3_Overlap_GENENAME.txt",
       na = '')

fwrite(genes_P13vsP7_D3_Untreated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D3/P13vsP7_D3_Untreated_only_ENSEMBL.txt",
       na = '')

fwrite(
  genes_P13vsP7_D3_Untreated_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP7_D3/P13vsP7_D3_Untreated_only_GENENAME.txt",
  na = ''
)

fwrite(genes_P13vsP7_D3_Treated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D3/P13vsP7_D3_Treated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_P13vsP7_D3_Treated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D3/P13vsP7_D3_Treated_only_GENENAME.txt",
       na = '')

fwrite(genes_P13vsP7_D3_Overlap$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D3/P13vsP7_D3_Overlap_ENSEMBL.txt",
       na = '')

fwrite(genes_P13vsP7_D3_Overlap$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D3/P13vsP7_D3_Overlap_GENENAME.txt",
       na = '')

fwrite(genes_P13vsP5_D3_Untreated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D3/P13vsP5_D3_Untreated_only_ENSEMBL.txt",
       na = '')

fwrite(
  genes_P13vsP5_D3_Untreated_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP5_D3/P13vsP5_D3_Untreated_only_GENENAME.txt",
  na = ''
)

fwrite(genes_P13vsP5_D3_Treated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D3/P13vsP5_D3_Treated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_P13vsP5_D3_Treated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D3/P13vsP5_D3_Treated_only_GENENAME.txt",
       na = '')

fwrite(genes_P13vsP5_D3_Overlap$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D3/P13vsP5_D3_Overlap_ENSEMBL.txt",
       na = '')

fwrite(genes_P13vsP5_D3_Overlap$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D3/P13vsP5_D3_Overlap_GENENAME.txt",
       na = '')

## Between passages at day 5, TvsUT

fwrite(genes_P7vsP5_D5_Untreated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D5/P7vsP5_D5_Untreated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_P7vsP5_D5_Untreated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D5/P7vsP5_D5_Untreated_only_GENENAME.txt",
       na = '')

fwrite(genes_P7vsP5_D5_Treated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D5/P7vsP5_D5_Treated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_P7vsP5_D5_Treated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D5/P7vsP5_D5_Treated_only_GENENAME.txt",
       na = '')

fwrite(genes_P7vsP5_D5_Overlap$GENEID %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D5/P7vsP5_D5_Overlap_ENSEMBL.txt",
       na = '')

fwrite(genes_P7vsP5_D5_Overlap$GENENAME %>% list(),
       file = "./output/list_genes/venn_P7vsP5_D5/P7vsP5_D5_Overlap_GENENAME.txt",
       na = '')

fwrite(genes_P13vsP7_D5_Untreated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D5/P13vsP7_D5_Untreated_only_ENSEMBL.txt",
       na = '')

fwrite(
  genes_P13vsP7_D5_Untreated_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP7_D5/P13vsP7_D5_Untreated_only_GENENAME.txt",
  na = ''
)

fwrite(genes_P13vsP7_D5_Treated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D5/P13vsP7_D5_Treated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_P13vsP7_D5_Treated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D5/P13vsP7_D5_Treated_only_GENENAME.txt",
       na = '')

fwrite(genes_P13vsP7_D5_Overlap$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D5/P13vsP7_D5_Overlap_ENSEMBL.txt",
       na = '')

fwrite(genes_P13vsP7_D5_Overlap$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP7_D5/P13vsP7_D5_Overlap_GENENAME.txt",
       na = '')

fwrite(genes_P13vsP5_D5_Untreated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D5/P13vsP5_D5_Untreated_only_ENSEMBL.txt",
       na = '')

fwrite(
  genes_P13vsP5_D5_Untreated_only$GENENAME %>% list(),
  file = "./output/list_genes/venn_P13vsP5_D5/P13vsP5_D5_Untreated_only_GENENAME.txt",
  na = ''
)

fwrite(genes_P13vsP5_D5_Treated_only$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D5/P13vsP5_D5_Treated_only_ENSEMBL.txt",
       na = '')

fwrite(genes_P13vsP5_D5_Treated_only$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D5/P13vsP5_D5_Treated_only_GENENAME.txt",
       na = '')

fwrite(genes_P13vsP5_D5_Overlap$GENEID %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D5/P13vsP5_D5_Overlap_ENSEMBL.txt",
       na = '')

fwrite(genes_P13vsP5_D5_Overlap$GENENAME %>% list(),
       file = "./output/list_genes/venn_P13vsP5_D5/P13vsP5_D5_Overlap_GENENAME.txt",
       na = '')
