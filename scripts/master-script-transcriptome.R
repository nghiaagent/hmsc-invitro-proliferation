### Load all packages using script under ./scripts/pkg_load.R
# If more packages are needed, please add them to pkg_load.R

source(file = "./scripts/pkg_load.R")

### Import transcript quantification estimates for downstream analysis
### I changed some folders around. Input data now in ./input and excluded from
### git repo inside .gitignore

tx2gene <- read.csv("2-Input/Homo_sapiens.GRCh38.91_tx2gene.csv")
head(tx2gene, 5)
folder <- c("2-Input/cell_quants")
salmon.dir <-
  as.matrix(read.csv(file = "2-Input/quant_filenames.csv", sep = ",", header =
                       F))
salmon.files <- file.path(folder, salmon.dir, "quant.sf")
names(salmon.files) <-
  as.matrix(read.csv(file = "2-Input/names.csv", sep = ",", header = F))
all(file.exists(salmon.files))
txi <-
  tximport(
    salmon.files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE,
    countsFromAbundance = "lengthScaledTPM"
  )  ### NOTE : 4545 transcripts missing ####
names(txi)
head(txi$counts, 3)
write.csv(txi, file = "4-Output/geneSMART_tximport_matrix.csv")

# --------------------------------------------------------------------------------------------------------------------

y <- DGEList(txi$counts)
dim(y)
csvfile <- file.path("2-Input/Cell_sample_table.csv")
sampleTable <- read.csv(csvfile, row.names = 1)
y$samples$names <- sampleTable$name
y$samples$filename <- sampleTable$filename
y$samples$Sample_ID <- sampleTable$Sample_ID
y$samples$cell_line <- sampleTable$cell_line
y$samples$Passage <- sampleTable$Passage
y$samples$Day <- sampleTable$Day
y$samples$Treatment <- sampleTable$Treatment
y$samples$day <- sampleTable$Day
y$samples$id <- sampleTable$Sample_ID
y$samples$run_date <- sampleTable$run_date
y$samples
geneid <- rownames(y)
genes <-
  select(
    Homo.sapiens,
    key = geneid,
    keytype = "ENSEMBL",
    columns = c("SYMBOL", "GENEID", "GENENAME", "TXCHROM"),
    multiVals = "first"
  )
geneid <- rownames(y)
ensembl91 <-
  useMart(host = "dec2017.archive.ensembl.org",
          biomart = "ENSEMBL_MART_ENSEMBL",
          dataset = "hsapiens_gene_ensembl")
attributes <- listAttributes(ensembl91)
attributes[1:20,]
genes <-
  select(
    ensembl91,
    keys = geneid,
    keytype = "ensembl_gene_id",
    columns = c(
      "ensembl_gene_id",
      "external_gene_name",
      "description",
      "entrezgene",
      "chromosome_name",
      "gene_biotype"
    )
  )
(colnames(genes) <-
    c("ENSEMBL", "SYMBOL", "GENENAME", "GENEID", "TXCHROM", "BIOTYPE"))
genes <- genes[!duplicated(genes$ENSEMBL),]
y$genes <- genes
head(genes, 5)
dim(genes)

# ------------------------------------------PLOT FOR SIZES-------------------------------------------------------------------------

png(
  "3-QC/Barplot of library sizes.png",
  width = 90,
  height = 30,
  units = 'cm',
  res = 300
)
col <- as.numeric(y$sample$timepoint)
dt <- colSums((y$counts) * 1e-6)
barplot(
  dt,
  names = colnames(dt),
  col = c("lightsalmon", "lightcoral", "steelblue1", "cornsilk")[col],
  las = 2,
  cex.names = 0.8
)
abline(h = 5, col = "black", lty = 3)
abline(h = 10, col = "black", lty = 3)
abline(h = 15, col = "black", lty = 3)
title(main = "Barplot of library sizes", ylab = "Library size (millions)")
dev.off()
# -----------------------------------------------------------------------------------------------------------------------

summary(dt)
CountMeans <- rowMeans(y$counts)
TotalCounts <- sum(CountMeans)
PercentReads <- (CountMeans / TotalCounts) * 100
AvgCounts <- rowMeans(y$counts)
Symbol <- y$genes$SYMBOL
GeneName <- y$genes$GENENAME
AvgCountsTable <-
  (data.frame(Symbol, GeneName, AvgCounts, PercentReads))
AvgCountsTable <- AvgCountsTable[order(-AvgCountsTable$AvgCounts),]
head(AvgCountsTable, 20)
write.csv(AvgCountsTable, file = "4-Output/Average_Counts.csv")
top5 <-
  ((sum(head(AvgCountsTable, 5)$AvgCounts)) / TotalCounts) * 100
top10 <-
  ((sum(head(AvgCountsTable, 10)$AvgCounts)) / TotalCounts) * 100
top25 <-
  ((sum(head(AvgCountsTable, 25)$AvgCounts)) / TotalCounts) * 100
dim(y)
summary(rowMeans(y$counts))
lcpm.AvgCounts <- cpm(AvgCounts, log = TRUE)
# ----------------------------------------------------PLOT-------------------------------------------------------------------------------
png(
  "3-QC/Read Count Distribution-B4.png",
  width = 45,
  height = 25,
  units = 'cm',
  res = 600
)
par(mar = c(5, 6, 4, 1) + .1)
hist(
  lcpm.AvgCounts,
  col = "salmon",
  border = "salmon",
  cex.lab = 2.5,
  cex.main = 3,
  cex.axis = 2,
  xlab = "Median log2-CPM",
  ylab = "No. of Transcripts",
  breaks = 100,
  xlim = c(-10, 20),
  main = "Read Count Distribution (before)"
)
dev.off()

# ------------------------------------------------------------------------------------------------------------------------------------

table(AvgCounts >= 1)
table(AvgCounts >= 10)
table(AvgCounts >= 100)
table(AvgCounts >= 1000)
table(AvgCounts >= 10000)
table(AvgCounts >= 100000)
table(rowSums(y$counts == 0) == 211)

median_cpm <- apply(cpm(y), 1, median)
expr_cutoff <- 0.5 # in cpm sum(median_cpm > expr_cutoff)
y.Filt <- y[median_cpm > expr_cutoff, ]
lcpm.Filt <- cpm(y.Filt$counts, log = TRUE)
log.cutoff <- log2(expr_cutoff)
summary(rowMeans(y.Filt$counts))
# -------------------------------------------------------------------------------------------------------
png(
  "3-QC/Density of count values.png",
  width = 10,
  height = 20,
  units = 'cm',
  res = 600
)
nsamples <- ncol(y)
col <- rainbow(nsamples)
par(mfrow = c(2, 1))
lcpm.Raw <- cpm(y$counts, log = TRUE)
plot(
  density(lcpm.Raw[, 1]),
  col = col[1],
  xlim = c(-10, 20),
  ylim = c(0, 0.3),
  main = "",
  xlab = ""
)
for (i in 2:nsamples) {
  den <- density(lcpm.Raw[, i])
  lines(den$x, den$y, col = col[i])
}
abline(
  v = log.cutoff,
  col = "red",
  lwd = 1,
  lty = 2,
  main = ""
)
title("Raw data", xlab = "log2-CPM")

lcpm.Filt <- cpm(y.Filt$counts, log = TRUE)
plot(
  density(lcpm.Filt[, 1]),
  col = col[1],
  xlim = c(-10, 20),
  ylim = c(0, 0.3),
  main = "",
  xlab = ""
)
for (i in 2:nsamples) {
  den <- density(lcpm.Filt[, i])
  lines(den$x, den$y, col = col[i])
}
abline(
  v = log.cutoff,
  col = "red",
  lwd = 1,
  lty = 2,
  main = ""
)
title("Filtered data (median CPM > 0.5)", xlab = "log2-CPM")
dev.off()
# ---------------------------------------------------------------------------------------------------------
AvgCounts <- rowMeans(y.Filt$counts)
lcpm.AvgCounts <- cpm(AvgCounts, log = TRUE)

png(
  "3-QC/Read Count Distribution-AFTER.png",
  width = 45,
  height = 25,
  units = 'cm',
  res = 600
)
par(mar = c(5, 6, 4, 1) + .1)
hist(
  lcpm.AvgCounts,
  col = "salmon",
  border = "salmon",
  cex.lab = 2.5,
  cex.main = 3,
  cex.axis = 2,
  xlab = "Median log2-CPM",
  ylab = "No. of Transcripts",
  breaks = 100,
  xlim = c(-10, 20),
  main = "Read Count Distribution (after)"
)
dev.off()

png(
  "3-QC/Sample heatmap.png",
  width = 60,
  height = 60,
  units = 'cm',
  res = 300
)
par(mfrow = c(1, 2))
lcpm.Raw <- cpm(y$counts, log = TRUE)
heatmap(cor(lcpm.Raw))
title("Raw data")
lcpm.Filt <- cpm(y.Filt$counts, log = TRUE)
heatmap(cor(lcpm.Filt))
title("Filtered data")
dev.off()
# ---------------------------------------------------------------------------------------------------------------
y.Norm <- calcNormFactors(y.Filt, method = "TMM")
y.Norm$samples$norm.factors
summary(y.Norm$samples$norm.factors)

png(
  "3-QC/BoxplotBefore.png",
  width = 90,
  height = 45,
  units = 'cm',
  res = 300
)
par(mfrow = c(2, 1))
lcpm.Filt <- cpm(y.Filt, log = TRUE)
col <- as.numeric(y$sample$tDay)
boxplot(
  lcpm.Filt,
  las = 2,
  col = c("lightsalmon", "lightcoral", "steelblue1", "cornsilk")[col],
  main = ""
)
abline(h = median(lcpm.Filt),
       col = "black",
       lty = 3)
title(main = "Unnormalised data", ylab = "log-counts")

png(
  "3-QC/BoxplotAfter.png",
  width = 90,
  height = 45,
  units = 'cm',
  res = 300
)
par(mfrow = c(2, 1))
lcpm.Norm <- cpm(y.Norm, log = TRUE)
boxplot(
  lcpm.Norm,
  las = 2,
  col = c("lightsalmon", "lightcoral", "steelblue1", "cornsilk")[col],
  main = ""
)
abline(h = median(lcpm.Norm),
       col = "black",
       lty = 3)
title(main = "TMM-Normalised data", ylab = "log-counts")
dev.off()
# -----------------------------------------------------THIS DOESNT WORK BUT I DONT THINK IT MATTERS-------------------------------------------------------------

dir.create("5-Glimma")
glMDSPlot(
  y.Norm,
  top = 500,
  labels = y.Norm$samples$subject,
  groups = y.Norm$samples[, c(1:9)],
  launch = TRUE,
  path = "5-Glimma",
  folder = "glMDSPlot",
  html = "MDS_plot"
)

png(
  "3-QC/MDS Plot_names.png",
  width = 15,
  height = 30,
  units = 'cm',
  res = 600
)
par(mfrow = c(2, 1))
plotMDS(
  y.Norm,
  top = 500,
  cex = 0.8,
  labels = y.Norm$samples$cell_line,
  col = as.numeric(y.Norm$samples$Treatment),
  main = "MDS Plot (SampleID)"
)
legend(
  "topright",
  legend = c("hMSC-20176", "hMSC-21558"),
  cex = 0.8,
  col = 1:16,
  pch = 16
)

plotMDS(
  y.Norm,
  top = 500,
  cex = 1,
  pch = 21,
  col = "white",
  bg = as.numeric(y.Norm$samples$Treatment),
  main = "MDS Plot (Symbol)"
) pch = x
0 = open square
1 = open circle
15 = solid square,
16 = solid circle, 21 = filled circle / border
legend(
  "topright",
  legend = c("hMSC-20176", "hMSC-21558"),
  cex = 0.8,
  col = 1:16,
  pch = 16
) dev.off()
png(
  "3-QC/MDPlots-1LI00.png",
  width = 11,
  height = 20,
  units = 'cm',
  res = 200
) par(mfrow = c(2, 1)) plotMD(y.Filt, column = 1, main = "First sample (raw)")
abline(
  h = 0,
  col = "red",
  lty = 2,
  lwd = 2
) plotMD(y.Norm, column = 1, main = "First sample (TMM-normalised)")
abline(
  h = 0,
  col = "red",
  lty = 2,
  lwd = 2
) dev.off()

# -------------------------------------------------------------------------------------------------------------------------------------------------------------

data <- y$samples
group <- paste(data$Treatment, data$Day, sep = ".")
group2 <- paste(data$cell_line, data$Treatment, sep = ".")
group3 <- paste(data$cell_line, data$Treatment, data$Day, sep = ".")
group4 <- paste(data$Treatment, data$Passage, sep = ".")
group5 <-
  paste(data$cell_line, data$Treatment, data$Day, data$Passage, sep = ".")

# -------------------------------- TREATED VS UNTREATED * POOLED PASSAGE * POOLED DAYS * POOLED POPULATION------------------------------------------------------
design <- model.matrix(~ 0 + data$Treatment)
colnames(design) <- c("Treated", "Untreated")

# ---------------------------------TREATED VS UNTREATED * POOLED POPULATIONS * POOLED PASSAGE * DIFFERENT DAYS---------------------------------------------------
design <- model.matrix(~ 0 + group)
colnames(design2) <-
  c("Treated_D3", "Treated_D5", "Untreated_D3", "Untreated_D5")
contrast.matrix <-
  makeContrasts("Treated_D3-Untreated_D3", levels = design2)

# ---------------------------------TREATED VS UNTREATED * SEPARATE POPULATIONS * POOLED PASSAGES * POOLED DAYS-------------------------------------------------
design <- model.matrix(~ 0 + group2)
colnames(design3) <-
  c(
    "hMSC 220176 Treated",
    "hMSC 220176 Untreated",
    "hMSC 221558 Treated",
    "hMSC 221558 Untreated"
  )

# -------------------------------_TREATED VS UNTREATED * SEP POPULATIONS * POOLED PASSAGES * SEP DAYS------------------------------------------------------
design <- model.matrix(~ 0 + group3)
colnames(design4) <-
  c(
    "hMSC 220176 Treated Day 3",
    "hMSC 220176 Treated Day 5",
    "hMSC 220176 Untreated Day 3",
    "hMSC 220176 Untreated Day 5",
    "hMSC 21558 Treated Day 3",
    "hMSC 21558 Treated Day 5",
    "hMSC 21558 Untreated Day 3",
    "hMSC 21558 Untreated Day 5"
  )

# -------------------------------TREATED VS UNTREATED * POOLED POP * SEP PASSAGES * POOLED DAYS------------------------------------------------------------

design <- model.matrix(~ 0 + group4)
colnames(design) <- c(
  "TreatedP13",
  "TreatedP5",
  "TreatedP7",
  "UntreatedP13",
  "UntreatedP5",
  "UntreatedP7"
)

# ---------------------------------EACH SAMPLE SEPARATE (DAYS, POP, PASSAGE, HEP NO HEP)-------------------------------------------------------------------

design6 <- model.matrix(~ 0 + group5)
#### Leaving names as it is... Theres so much....

# ------------------------------------------------------------------------------------------------------------------------------------------------------


# -----------------------------------------------FOLLOWING NICK-----------------------------------------------------------------------
png(
  "3-QC/SA Plot.png",
  width = 30,
  height = 15,
  units = 'cm',
  res = 600
)
par(mfrow = c(1, 2))
v <- voom(y.Norm, design, plot = TRUE)
plotSA(fit2, main = "Final model: Mean−variance trend")
dev.off()

v2 <- voom(y.Norm, design)
corfit <-
  duplicateCorrelation(v2, design, block = y.Norm$samples$cell_line)
v <-
  voom(
    y.Norm,
    design,
    block = y.Norm$samples$cell_line,
    correlation =
      corfit$consensus
  )
fit2 <-
  lmFit(v,
        design,
        block = y.Norm$samples$cell_line,
        correlation = corfit$consensus)
contrast.matrix1 <-
  makeContrasts("TreatedP5-TreatedP13", levels = design)
contrast.matrix2 <-
  makeContrasts("UntreatedP5-UntreatedP13", levels = design)
fit2c <- contrasts.fit(fit2, contrast.matrix1)
fit3c <- contrasts.fit(fit2, contrast.matrix2)
names(fit2c)
names(fit3c)

# METHOD1
fit2c <- eBayes(fit2c)
# Then use TopTable
# ---------------------------------------------------------------------------------------------------
# METHOD2
lfc <- log2(1.1)
fit2c <- treat(fit2c, lfc = lfc)
fit3c <- treat(fit3c, lfc = lfc)
save(fit2c, file = "4-Output/v.rda")
save(fit3c, file = "4-Output/v2.rda")
# Then use topTreat

results1 <- decideTests(fit2c)
results2 <- decideTests(fit3c)
summary(results1)
summary(results2)
vennCounts(results1)
vennCounts(results2)
ngenes1 <- which(results1[, 1] != 0)
ngenes1 <- length(ngenes1)
ngenes2 <- which(results2[, 1] != 0)
ngenes2 <- length(ngenes2)

png(
  "4-Output/Venn DiagramP5-P13+Hep.png",
  width = 30,
  height = 15,
  units = 'cm',
  res = 600
)
par(mfrow = c(1, 2))
vennDiagram(
  results1,
  include = c("both"),
  circle.col = c("blue", "yellow", "green", "orange"),
  counts.col = c("blue3"),
  cex = c(1, 0.8, 0.8)
)
vennDiagram(
  results1,
  include = c("up", "down"),
  circle.col = c("blue", "yellow", "green", "orange"),
  counts.col = c("red", "green3"),
  cex = c(1, 0.7, 0.7)
)
mtext(
  "hMSC",
  side = 3,
  line = -2,
  outer = TRUE,
  cex = 1.8
)
mtext(
  paste0("(", ngenes1, " genes)"),
  side = 3,
  line = -3.5,
  outer = TRUE,
  cex = 1.5
)
dev.off()

png(
  "4-Output/Venn DiagramP5-P13UT.png",
  width = 30,
  height = 15,
  units = 'cm',
  res = 600
)
par(mfrow = c(1, 2))
vennDiagram(
  results2,
  include = c("both"),
  circle.col = c("blue", "yellow", "green", "orange"),
  counts.col = c("blue3"),
  cex = c(1, 0.8, 0.8)
)
vennDiagram(
  results2,
  include = c("up", "down"),
  circle.col = c("blue", "yellow", "green", "orange"),
  counts.col = c("red", "green3"),
  cex = c(1, 0.7, 0.7)
)
mtext(
  "hMSC",
  side = 3,
  line = -2,
  outer = TRUE,
  cex = 1.8
)
mtext(
  paste0("(", ngenes2, " genes)"),
  side = 3,
  line = -3.5,
  outer = TRUE,
  cex = 1.5
)
dev.off()
# -----------------------------------------------FOLLOWING ANOTHER PROTOCOL IF NEEDED ----------------------------------------------------------
design <- model.matrix(~ group4)
dgeObj <- estimateCommonDisp(y.Norm, verbose = TRUE)
# Disp = 0.30772 , BCV = 0.5547
dgeObj <- estimateGLMTrendedDisp(dgeObj)
dgeObj <- estimateTagwiseDisp(dgeObj)
plotBCV(
  dgeObj,
  xlab = "Average log CPM",
  ylab = "Biological coefficient of variation",
  pch = 16,
  cex = 0.2,
  col.common = "red",
  col.trend = "blue",
  col.tagwise = "black"
)
fitd <- glmFit(dgeObj, design)
names(fitd)
head(coef(fitd))
lrt.BvsL <- glmLRT(fitd)
topTags(lrt.BvsL)

PvsV <-
  makeContrasts(group4Treated.P5 - group4Untreated.P5, levels = design)
lrt.pVsV <- glmLRT(fitd, contrast = PvsV)
topTags(lrt.pVsV)
results <- as.data.frame(topTags(lrt.pVsV, n = Inf))
dim(results)
summary(de <- decideTestsDGE(lrt.pVsV))
detags <- rownames(dgeObj)[as.logical(de)]
plotSmear(lrt.pVsV, de.tags = detags)
# ------------------------------------------------------------------------------------------------------------------------------------------------------------

png(
  "4-Output/MDplot-Average.png",
  width = 20,
  height = 30,
  units = 'cm',
  res = 600
)
plotMD(
  fit2c,
  coef = 1,
  status = results[, 1],
  values = c(-1, 1),
  main = colnames(fit2c)[1]
)
dev.off()

glMDPlot(
  fit2c,
  coef = 1,
  status = results[, 1],
  main = colnames(fit2c)[1],
  side.main = "SYMBOL",
  counts = y.Norm$counts,
  groups = v$targets$chip,
  launch = TRUE,
  path = "5-Glimma",
  folder = "glMDPlots",
  html = "Treated P5 - Treated P13"
)
glMDPlot(
  fit3c,
  coef = 1,
  status = results[, 1],
  main = colnames(fit3c)[2],
  side.main = "SYMBOL",
  counts = y.Norm$counts,
  groups = v$targets$chip,
  launch = TRUE,
  path = "5-Glimma",
  folder = "glMDPlots",
  html = "Untreated P5 - Untreated P13"
)


### This shows an error but seems to be working####

# --------------------------------------------------CAN ALSO PLOT MORE THAN ONE COMPARISON BUT ADD PAR--------------------------------------------------------

png(
  "4-Output/Volcano_plots.png",
  width = 20,
  height = 30,
  units = 'cm',
  res = 600
)
cutoff = -log10(0.05)
volcanoplot(
  fit2c,
  coef = 1,
  highlight = 25,
  names = v$genes$SYMBOL,
  main = colnames(fit2c)[1]
)
dev.off()

png(
  "4-Output/Volcano_plots2.png",
  width = 20,
  height = 30,
  units = 'cm',
  res = 600
)
cutoff = -log10(0.05)
volcanoplot(
  fit3c,
  coef = 1,
  highlight = 25,
  names = v$genes$SYMBOL,
  main = colnames(fit2c)[1]
)
dev.off()

glXYPlot(
  x = fit2c$coef[, 1],
  y = fit2c$p.value[, 1],
  xlab = "logFC",
  ylab = "p-value",
  status = results[, 1],
  main = colnames(fit2c)[1],
  side.main = "SYMBOL",
  anno = fit2c$genes,
  counts = y.Norm$counts,
  groups = v$targets$timepoint,
  path = "5-Glimma",
  folder = "glVolcano",
  html = "Treated P5 - Treated P13"
)



# ---------------------------------------------------ADD MORE-----------------------------------------------------------------------------------------------
par(mfrow = c(3, 1))
volcanoplot(
  fit2,
  coef = 2,
  highlight = 25,
  names = v$genes$SYMBOL,
  main = colnames(fit2)[2]
)
volcanoplot(
  fit2,
  coef = 3,
  highlight = 25,
  names = v$genes$SYMBOL,
  main = colnames(fit2)[3]
)
# ------------------------------------------------------------CAN ALSO START FROM HERE NEXT TIME----------------------------------------------------------------------------------------------


r1 <- topTreat(fit2c,
               adjust = "BH",
               coef = 1,
               n = Inf)
write.csv(r1, file = "4-Output/topTable_Treated.csv")
sig <- r1$adj.P.Val < 0.05
cat("No.Sig.Genes.Treated:", length(which(sig == 1)))

r2 <- topTreat(fit3c,
               adjust = "BH",
               coef = 1,
               n = Inf)
write.csv(r2, file = "4-Output/topTable_Untreated.csv")
sig2 <- r2$adj.P.Val < 0.05
cat("No.Sig.Genes.Untreated:", length(which(sig2 == 1)))
# ---------------------------------------------------GENE SET TESTING-------------------------------------------------------------------------------------------------------

results <- as.data.frame(r1)
results2 <- as.data.frame(r2)
genes <- as.list(r1$adj.P.Val < 0.05)
names(genes) <- results$SYMBOL
result <- as.data.frame(genes[1:265])
result <- names(result)
gene.df <- bitr(
  result,
  fromType = "SYMBOL",
  toType = c("ENSEMBL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)
gene_symbols <- unlist(gene.df[, "ENTREZID"])
ggo <- groupGO(
  gene     = gene_symbols,
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  level    = 3,
  readable = TRUE
)
ggo_df <- data.frame(ggo)
write.csv(ggo_df, file = "4-Output/GO_Analysis.csv")

genes <- as.list(r2$adj.P.Val < 0.05)
names(genes) <- results2$SYMBOL
result2 <- as.data.frame(genes[1:265])
result2 <- names(result)
gene.df <- bitr(
  result2,
  fromType = "SYMBOL",
  toType = c("ENSEMBL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)
gene_symbols <- unlist(gene.df[, "ENTREZID"])
ggo <- groupGO(
  gene     = gene_symbols,
  OrgDb    = org.Hs.eg.db,
  ont      = "BP",
  level    = 3,
  readable = TRUE
)
ggo_df <- data.frame(ggo)
write.csv(ggo_df, file = "4-Output/GO_Analysis2.csv")

# ----------------------------------------------------------------
results.ord <- r1[order(-r1[, "logFC"]), ]
head(results.ord)
ranks <- results.ord$logFC
names(ranks) <- results.ord$SYMBOL
head(ranks)

results.ord2 <- r2[order(-r2[, "logFC"]), ]
head(results.ord2)
ranks2 <- results.ord2$logFC
names(ranks2) <- results.ord2$SYMBOL
head(ranks2)

pathways.hallmark <- gmtPathways("2-Input/h.all.v7.0.symbols.gmt")
head(pathways.hallmark)
fgseaRes <-
  fgseaMultilevel(pathways = pathways.hallmark, stats = ranks)
fgseaRes2 <-
  fgseaMultilevel(pathways = pathways.hallmark, stats = ranks2)

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
gene.in.pathway <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  inner_join(r1, by = "SYMBOL")

fgseaResTidy2 <- fgseaRes2 %>%
  as_tibble() %>%
  arrange(desc(NES))
gene.in.pathway2 <- pathways.hallmark %>%
  enframe("pathway", "SYMBOL") %>%
  unnest(cols = c(SYMBOL)) %>%
  inner_join(r1, by = "SYMBOL")


dev.new()
fgseaResTidy$adjPvalue <-
  ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey",
          "significant" = "red")
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = factor(adjPvalue))) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Hallmark pathways Enrichment Score from GSEA")
dev.off()
dev.new()
fgseaResTidy2$adjPvalue <-
  ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey",
          "significant" = "red")
ggplot(fgseaResTidy2, aes(reorder(pathway, NES), NES, fill = factor(adjPvalue))) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  coord_flip() +
  labs(x = "Pathway", y = "Normalized Enrichment Score",
       title = "Hallmark pathways Enrichment Score from GSEA")
dev.off()


plotEnrichment(pathway = pathways.hallmark[["HALLMARK_HEDGEHOG_SIGNALING"]], ranks)

plotGseaTable(pathways.hallmark[fgseaRes$pathway[fgseaRes$padj < 0.05]], ranks, fgseaRes,
              gseaParam = 0.5)

plotEnrichment(pathway = pathways.hallmark[["HALLMARK_HEDGEHOG_SIGNALING"]], ranks2)

plotGseaTable(pathways.hallmark[fgseaRes2$pathway[fgseaRes2$padj < 0.05]], ranks2, fgseaRes,
              gseaParam = 0.5)


sig.path <-
  fgseaResTidy$pathway[fgseaResTidy$adjPvalue == "significant"]
sig.gen <-
  unique(na.omit(gene.in.pathway$SYMBOL[gene.in.pathway$pathway %in% sig.path]))
h.dat <-
  dcast(gene.in.pathway[, c(1, 2)], SYMBOL ~ pathway, drop = FALSE)
symbols <- h.dat$SYMBOL
rownames(h.dat) <- symbols
h.dat <- h.dat[, -1, drop = FALSE]
rownames(h.dat)
rownames(h.dat) <- symbols
rownames(h.dat)
h.dat <- h.dat[rownames(h.dat) %in% sig.gen, , drop = FALSE]
h.dat <- cbind(row.names(h.dat), h.dat)
h.dat <- h.dat[, colnames(h.dat) %in% sig.path]
table(data.frame(rowSums(h.dat)))
h.dat <- h.dat[data.frame(rowSums(h.dat)) >= 3, ]
topTable <- r1[r1$SYMBOL %in% rownames(h.dat), ]

topTableAligned <-
  topTable[which(topTable$SYMBOL %in% rownames(h.dat)),]
topTableAligned <-
  topTableAligned[match(rownames(h.dat), topTableAligned$SYMBOL),]
all(topTableAligned$SYMBOL == rownames(h.dat))

dfMinusLog10FDRGenes <-
  data.frame(-log10(topTableAligned[which(topTableAligned$SYMBOL %in% rownames(h.dat)), 'adj.P.Val']))
dfMinusLog10FDRGenes[dfMinusLog10FDRGenes == 'Inf'] <- 0

dfFoldChangeGenes <-
  data.frame(topTableAligned[which(topTableAligned$SYMBOL %in% rownames(h.dat)), 'logFC'])
dfGeneAnno <- data.frame(dfMinusLog10FDRGenes, dfFoldChangeGenes)
colnames(dfGeneAnno) <- c('Gene score', 'LogFC')
dfGeneAnno[, 2] <- ifelse(
  dfGeneAnno$LogFC > 0,
  'Up-regulated',
  ifelse(dfGeneAnno$LogFC < 0, 'Down-regulated', 'Unchanged')
)

colours <- list('LogFC' = c(
  'Up-regulated' = 'royalblue',
  'Down-regulated' = 'yellow'
))
haGenes <- rowAnnotation(
  df = dfGeneAnno,
  col = colours,
  width = unit(1, 'cm'),
  annotation_name_side = 'top'
)

dfEnrichment <- fgseaRes[, c("pathway", "NES")]
dfEnrichment <-
  dfEnrichment[dfEnrichment$pathway %in% colnames(h.dat)]
dd <- dfEnrichment$pathway
dfEnrichment <- dfEnrichment[, -1]
rownames(dfEnrichment) <- dd
colnames(dfEnrichment) <- 'Normalized\n Enrichment score'
haTerms <- HeatmapAnnotation(
  df = dfEnrichment,
  Term = anno_text(
    colnames(h.dat),
    rot = 45,
    just = 'right',
    gp = gpar(fontsize = 12)
  ),
  annotation_height = unit.c(unit(1, 'cm'), unit(8, 'cm')),
  annotation_name_side = 'left'
)

hmapGSEA <- Heatmap(
  h.dat,
  name = 'GSEA hallmark pathways enrichment',
  split = dfGeneAnno[, 2],
  col = c('0' = 'white', '1' = 'forestgreen'),
  rect_gp = gpar(col = 'grey85'),
  cluster_rows = TRUE,
  show_row_dend = TRUE,
  row_title = 'Top Genes',
  row_title_side = 'left',
  row_title_gp = gpar(fontsize = 11, fontface = 'bold'),
  row_title_rot = 90,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 11, fontface = 'bold'),
  row_names_side = 'left',
  row_dend_width = unit(35, 'mm'),
  cluster_columns = TRUE,
  show_column_dend = TRUE,
  column_title = 'Enriched terms',
  column_title_side = 'top',
  
  column_title_gp = gpar(fontsize = 12, fontface = 'bold'),
  column_title_rot = 0,
  show_column_names = FALSE,
  show_heatmap_legend = FALSE,
  clustering_distance_columns = 'euclidean',
  clustering_method_columns = 'ward.D2',
  clustering_distance_rows = 'euclidean',
  clustering_method_rows = 'ward.D2',
  bottom_annotation = haTerms
)
dev.on()

tiff(
  "GSEA_enrichment_2.tiff",
  units = "in",
  width = 13,
  height = 22,
  res = 400
)
draw(
  hmapGSEA + haGenes,
  heatmap_legend_side = 'right',
  annotation_legend_side = 'right'
)
dev.off()
