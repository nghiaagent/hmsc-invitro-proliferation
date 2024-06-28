# Perform DGE beforehand
# Load 6_post_camera_function
# This performs camera on logCPM (quant_DGE_voom)
# Automatically exports results

enrichments_camera <- purrr::map(1:ncol(matrix_contrasts), \(x) run_camera(fit_contrasts, x))

names(enrichments_camera) <- colnames(matrix_contrasts)

saveRDS(enrichments_camera, file.path('output',
                                      'data_enrichment',
                                      'camera',
                                      "camera_all.RDS"))

# Barcode plots of selected pathways

barcodeplot(topTable(fit_contrasts,
                     coef = 1,
                     number = Inf,
                     sort.by = "none")$t,
            index = msigdb_reactome$REACTOME_OLFACTORY_SIGNALING_PATHWAY)

glimmaVolcano(x = fit_contrasts[msigdb_reactome$REACTOME_OLFACTORY_SIGNALING_PATHWAY,],
              counts = quant_DGE_voom[msigdb_reactome$REACTOME_OLFACTORY_SIGNALING_PATHWAY,]$E,
              groups = quant_DGE_voom[msigdb_reactome$REACTOME_OLFACTORY_SIGNALING_PATHWAY,]$targets$condition_ID,
              xlab = "logFC",
              transform.counts = "none",
              coef = 1)