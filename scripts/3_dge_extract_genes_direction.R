# Extract genes based on their direction of regulation between P5 - P7 and P7 - P13
# (e.g. up - up, down - down, down - up, up - down)

tests <- decideTests(fit_contrasts)

genes_up_up <-
  quant_DGE_voom$genes[which(tests[, 13] == 1 & tests[, 14] == 1),]

genes_up_down <-
  quant_DGE_voom$genes[which(tests[, 13] == 1 & tests[, 14] == -1),]

genes_down_up <-
  quant_DGE_voom$genes[which(tests[, 13] == -1 & tests[, 14] == 1),]

genes_down_down <-
  quant_DGE_voom$genes[which(tests[, 13] == -1 & tests[, 14] == -1),]

# Extract genes based on Hep response 

# Export lists

fwrite(
  genes_up_up$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_up_ENSEMBL.txt",
  na = ''
)

fwrite(
  genes_up_up$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_up_GENENAME.txt",
  na = ''
)


fwrite(
  genes_up_down$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_down_ENSEMBL.txt",
  na = ''
)

fwrite(
  genes_up_down$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_down_GENENAME.txt",
  na = ''
)

fwrite(
  genes_up_down$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_down_ENSEMBL.txt",
  na = ''
)

fwrite(
  genes_up_down$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_up_down_GENENAME.txt",
  na = ''
)

fwrite(
  genes_down_down$GENEID %>% list(),
  file = "./output/list_genes/betweenpassages/genes_down_down_ENSEMBL.txt",
  na = ''
)

fwrite(
  genes_down_down$GENENAME %>% list(),
  file = "./output/list_genes/betweenpassages/genes_down_down_GENENAME.txt",
  na = ''
)