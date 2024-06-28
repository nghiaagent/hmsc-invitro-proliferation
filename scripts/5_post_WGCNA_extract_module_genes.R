# Load data

## GCN with all samples

load(file.path(
  '.',
  'output',
  'data_WGCNA',
  'WGCNA_allsamples',
  'GCN_output.Rdata'
))

## Add module assignment to gene list

GCN$genes$module <- GCN$net$colors %>%
  labels2colors()

GCN$genes <- GCN$genes %>%
  mutate(ENTREZID = as.character(ENTREZID))

## Extract gene lists
## ENTREZ for GSVA
## GENENAME for String/PPI

GCN_genelists_ENTREZ <- split(GCN$genes, GCN$genes$module, drop = TRUE) %>%
  lapply(function(x) {
    select(x, ENTREZID) %>%
      as_vector %>%
      unname()
  })

GCN_genelists_GENENAME <- split(GCN$genes, GCN$genes$module, drop = TRUE) %>%
  lapply(function(x) {
    select(x, GENENAME) %>%
      as_vector %>%
      unname()
  })

names(GCN_genelists_ENTREZ)   <- str_c("WGCNA_allsamples_", names(GCN_genelists_ENTREZ))
names(GCN_genelists_GENENAME) <- str_c("WGCNA_allsamples_", names(GCN_genelists_GENENAME))

# Define function to write list to GMT file
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

# Write gene list in ENTREZID to GMT file

writeGMT(
  GCN_genelists_ENTREZ,
  file.path(
    ".",
    "input",
    "genesets",
    "GCN_sets_WGCNA_allsamples_only.gmt"
  )
)

# Write gene list with GENENAME to text files for STRING

fwrite(
  GCN_genelists_GENENAME[["WGCNA_allsamples_turquoise"]] %>% list(),
  file = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "turquoise_genes.txt"
  ),
  na = ""
)

fwrite(
  GCN_genelists_GENENAME[["WGCNA_allsamples_red"]] %>% list(),
  file = file.path(
    ".",
    "output",
    "data_WGCNA",
    "WGCNA_allsamples",
    "red_genes.txt"
  ),
  na = ""
)