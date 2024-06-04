# Load data

## GCN with all samples

load(file.path(
  '.',
  'output',
  'data_WGCNA',
  'WGCNA_allsamples',
  'GCN_output.Rdata'
))

GCN_allsamples <- GCN

## GCN with only D3 UT samples

load(file.path('.',
               'output',
               'data_WGCNA',
               'WGCNA_D3_UT',
               'GCN_output.Rdata'))

GCN_subset <- GCN

rm(GCN)

## Add module assignment to gene list

GCN_subset$genes$module <- GCN_subset$net$colors %>%
  labels2colors()

GCN_subset$genes <- GCN_subset$genes %>%
  mutate(ENTREZID = as.character(ENTREZID))

GCN_allsamples$genes$module <- GCN_allsamples$net$colors %>%
  labels2colors()

GCN_allsamples$genes <- GCN_allsamples$genes %>%
  mutate(ENTREZID = as.character(ENTREZID))

## Extract gene lists

### subset

GCN_subset_genelists <- split(GCN_subset$genes,
                              GCN_subset$genes$module,
                              drop = TRUE) %>%
  lapply(function(x) {
    select(x, ENTREZID) %>%
      as_vector %>%
      unname()
  })

### All data

GCN_allsamples_genelists <- split(GCN_allsamples$genes,
                                  GCN_allsamples$genes$module,
                                  drop = TRUE) %>%
  lapply(function(x) {
    select(x, ENTREZID) %>%
      as_vector %>%
      unname()
  })

GCN_genelists_WGCNA <-
  c(GCN_allsamples_genelists, GCN_subset_genelists)

names(GCN_genelists_WGCNA) <- c(str_c("WGCNA_allsamples_",
                                      names(GCN_allsamples_genelists)),
                                str_c("WGCNA_subset_",
                                      names(GCN_subset_genelists)))

## Combine with CEMiTool files to make GCN gene lists

GCN_genelists_CEMi <- readGMT(file.path(".",
                                        "output",
                                        "data_WGCNA",
                                        "CEMiTool",
                                        "modules_genes.gmt"))

names(GCN_genelists_CEMi) <- c(str_c("CEMiTool_",
                                     names(GCN_genelists_CEMi)))

GCN_genelists <- c(GCN_genelists_WGCNA, GCN_genelists_CEMi)
  
  # Define function to write list to GMT file
  #Create a gmt (gene matrix transposed) file
  ### Createss a gmt (gene matrix transposed) file such as those
  ### provided by mSigDB or geneSigDB, from an R list object.
  ### Function by Levi Waldron.
  ### R list object that will be converted to GMT file.  Each element
  ### should contain a vector of gene names, and the names of the
  ### elements will used for the gene set names
  ### Output file name for .gmt file
  
writeGMT <- function (object, fname) {
  if (class(object) != "list")
    stop("object should be of class 'list'")
  
  if (file.exists(fname))
    unlink(fname)
  
  for (iElement in 1:length(object)) {
    write.table(
      t(c(make.names(rep(
        names(object)[iElement], 2
      )), object[[iElement]])),
      sep = "\t",
      quote = FALSE,
      file = fname,
      append = TRUE,
      col.names = FALSE,
      row.names = FALSE
    )
  }
}

writeGMT(GCN_genelists,
         file.path(".",
                   "input",
                   "genesets",
                   "GCN_sets.gmt"))

saveRDS(GCN_genelists,
        file = file.path('.',
                         'output',
                         'data_WGCNA',
                         'GCN_genelists.RDS'))