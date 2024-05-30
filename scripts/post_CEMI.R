
quant_DGE_CEMI <- quant_DGE_voom[, quant_DGE_voom$targets$Day == "D3" &
                                   quant_DGE_voom$targets$Treatment == "Untreated"]

temp <- quant_DGE_CEMI

rownames(temp$genes) <- temp$genes$GENEID

order <- order(fit$Amean, decreasing = TRUE)

temp <- temp[order, ]

temp <- temp[complete.cases(dplyr::select(temp$genes,
                            !msigdb_h)), ]

temp <- temp[!duplicated(temp$genes$ENTREZID), ]

rownames(temp$E) <- temp$genes$ENTREZID

rownames(temp$genes) <- temp$genes$ENTREZID

temp <- temp[order(temp$genes$GENEID, decreasing = FALSE), ]

temp$genes <- temp$genes %>% dplyr::select(!GENEID)

quant_DGE_CEMI <- temp

quant_DGE_CEMI$targets <- quant_DGE_CEMI$targets %>%
  dplyr::select(ID, Passage) %>%
  rename("SampleName" = "ID",
         "Class" = "Passage") %>%
  mutate(SampleName = as.character(SampleName))

cem <- cemitool(
  expr = quant_DGE_CEMI$E %>% as.data.frame(),
  annot = quant_DGE_CEMI$targets,
  gmt = read_gmt("./input/genesets/msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.entrez.gmt"),
  verbose = TRUE
)

save(cem, file = "./output/data_WGCNA/CEMiTool/CEMI_output.RData")

write_files(cem, "./output/data_WGCNA/CEMiTool/", force = TRUE)
