### Don't source this file by itself; call in from another file after running env_prep.R
### Source 
### Import salmon transcript quantification results to a single matrix, export that matrix.

# tx2gene is a df linking transcript ID to gene ID Hence, "transcript to gene"
# We use this to turn transcript quantification results to gene quant results.
# Import tx2gene df

tx2gene <- read.csv("./input/Homo_sapiens.GRCh38.91_tx2gene.csv")
head(tx2gene, 5)

# Construct list containing salmon transcript quantification files for each sample
# folder_tx_quant: root folder containing all quantification results folders
# dir_tx_quant: list of quantification results folders paths
# files_tx_quant: list of quantification results files paths

folder_tx_quant <- c("./input/cell_quants")

dir_tx_quant <-
  as.matrix(read.csv(file = "./Input/quant_filenames.csv",
                     sep = ",",
                     header = F))

files_tx_quant <-
  file.path(folder_tx_quant, 
            dir_tx_quant, 
            "quant.sf")

names(files_tx_quant) <-
  as.matrix(read.csv(file = "./input/names.csv", 
                     sep = ",", 
                     header = F))

all(file.exists(files_tx_quant)) # Make sure that all transcript files exist

# Import quantification results

matrix_tx_quant <-
  tximport(
    files_tx_quant,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE,
    countsFromAbundance = "lengthScaledTPM"
  )  ### NOTE : 4545 transcripts missing ####

names(matrix_tx_quant) # Check format

# Export quantification results

write.csv(matrix_tx_quant, file = "./output/matrix_transcript_quantification.csv")