camera_all <- purrr::map(
    1:ncol(contrasts),
    \(x) {
        run_camera(  
            elist = rlog_camera,
            design = design,
            contrasts = contrasts,
            genesets = list_gmt,
            coef = x,
            inter.gene.cor = 0.01,
            sort = TRUE
        )
    }
)

names(camera_all) <- colnames(contrasts)

saveRDS(
    enrichments_camera,
    here::here(
        "output",
        "data_enrichment",
        "camera",
        "camera_all.RDS"
    )
)                     

# 
# # Barcode plots of selected pathways
# 
# barcodeplot(topTable(fit_contrasts,
#                      coef = 15,
#                      number = Inf,
#                      sort.by = "none")$t,
#             index = msigdb_h$HALLMARK_TNFA_SIGNALING_VIA_NFKB)
# 
# glimmaVolcano(x = fit_contrasts[msigdb_reactome$REACTOME_OLFACTORY_SIGNALING_PATHWAY,],
#               counts = quant_DGE_voom[msigdb_reactome$REACTOME_OLFACTORY_SIGNALING_PATHWAY,]$E,
#               groups = quant_DGE_voom[msigdb_reactome$REACTOME_OLFACTORY_SIGNALING_PATHWAY,]$targets$condition_ID,
#               xlab = "logFC",
#               transform.counts = "none",
#               coef = 1)