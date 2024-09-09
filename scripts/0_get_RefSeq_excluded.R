# Get list of RefSeq IDs removed from RefSeq between design of AmpliSeq panel and now

ampliseq_manifest <- read_csv("./input/annotation/hg19_AmpliSeq_Transcriptome_21K_v1_DataSheet.csv",
                            skip = 8) %>%
  select(c("Name",
            "Gene_1",
            "Amplicon_ID"))

refseq_names <- read_tsv("./input/annotation/GCF_000001405.39_GRCh38.p13_rna_filter.fna.fai",
                         col_names = FALSE)

diff <- setdiff(ampliseq_manifest$Name, refseq_names$X1) %>%
  AnnotationDbi::select(
    org.Hs.eg.db,
    keys = .,
    keytype = "REFSEQ",
    columns = c("REFSEQ", "GENENAME", "SYMBOL", "ENTREZID", "GENETYPE")
  ) %>%
  arrange(SYMBOL)
