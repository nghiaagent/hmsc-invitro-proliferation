# This performs camera on logCPM (quant_DGE_voom)
# Automatically export results

enrichments_camera <- purrr::map(1:ncol(fit_contrasts$contrasts), \(x) run_camera(fit_contrasts, x, sort = TRUE))

names(enrichments_camera) <- colnames(fit_contrasts$contrasts)

saveRDS(enrichments_camera, file.path('output',
                                      'data_enrichment',
                                      'camera',
                                      "camera_all.RDS"))

