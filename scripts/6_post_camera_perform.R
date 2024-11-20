# Perform camera on all contrasts
# Load function script beforehand

camera_all <- purrr::map(
    seq_len(ncol(contrasts)),
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

# Save data

names(camera_all) <- colnames(contrasts)

saveRDS(
    camera_all,
    here::here(
        "output",
        "data_enrichment",
        "camera",
        "camera_all.RDS"
    )
)
