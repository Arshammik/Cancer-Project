# WGCNA Source Code
#==============================
#Chapter One
#DEG Evaluation

setwd("~/desktop/THCA")
ex <- read.delim("~/Desktop/THCA/ex_LUAD.tsv")
cn <- ex #backup

#cn3 <- aggregate(.~gene_name,cn,mean)
#can't run on local system for this scale, maybe on server

library(dplyr)
dim(cn)
cn <- distinct(cn, gene_id, .keep_all = TRUE)
dim(cn)
# lose 1233 genes but with geneID we loose nothing

rownames(cn) <- cn$gene_id
cn <- cn[,-1:-2]

# no idea why!!!
cn <- cn[,-c(178,223,226,246,339,428,485,474)]

#sample types
sampleID <- read.delim("../THCA/TCGA-LUAD/Data/gdc_sample_sheet.2022-06-11.tsv")
i1 <- sampleID$Data.Type == "Gene Expression Quantification"
sampleID <- sampleID[i1, , drop = F]
sampleID <- sampleID[,7:8]
sampleID <- distinct(sampleID, Sample.ID, .keep_all = T)
a <- colnames(cn)
for(i in 1:3){a <- sub('\\.', '-', a)}
b <- data.frame(IDs = a)
c <- merge(b, sampleID, by.x = "IDs", by.y = "Sample.ID", sort = F)
c$Sample.Type <- paste0(c$Sample.Type,".", seq_along(c$Sample.Type))
colnames (cn) <- c$Sample.Type

# we reached to our tidy data set

library(DESeq2)
gr <- factor(c$Sample.Type)
colSums(cn) #Through this we can demonstrates 
colData <- data.frame(group = gr, type = "paired-end")
cds <- DESeqDataSetFromMatrix(cn, colData, design = ~group)
cds <- DESeq(cds)
cnt <- log2(1+(counts(cds, normalize = T))) #getting normalized counts

#write.table(cnt, "~/desktop/LUAD/TCGA-LUAD/Results/expression(log2+1)(cnt).csv",  quote = F, col.names = T, row.names = T, sep = "\t")

#DEGs
dif <- results(cds, c("group", "Solid Tissue Normal" , "Primary Tumor"))
#write.table(dif, "~/desktop/Systematic Review/TCGA-LUAD/Results/dif.csv",  quote = F, col.names = T, row.names = T, sep = "\t")

#checkpoint
sorted_dif <- data.frame(results(cds, c("group", "Solid Tissue Normal" , "Primary Tumor"))) #which one should be the first — order is important
sorted_dif$padj <- p.adjust(sorted_dif$pvalue, method = "BH")
sorted_dif <- sorted_dif[order(sorted_dif$padj),]



setwd("~/desktop")
X  <- subset(sorted_dif, log2FoldChange > 1  & padj < 0.05)
Y  <- subset(sorted_dif, log2FoldChange < -1  & padj < 0.05)
write.table(subset(X), "Up regulated (LUAD-TCGA).csv",  quote = F, col.names = T, row.names = T, sep = "\t")
write.table(subset(Y), "Down regulated (LUAD-TCGA).csv",  quote = F, col.names = T, row.names = T, sep = "\t")

#==============================
#Chapter One
#WGCNA

data <- ex
dim(data)
data <- distinct(data, gene_name, .keep_all = TRUE)
dim(data)
rownames(data) <- data$gene_id
data <- data [,-1:-2]

# quality control
#################
library(WGCNA)
gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

# trim
data <- data[gsg$goodGenes == TRUE,]
gsg <- goodSamplesGenes(t(data))
gsg$allOK    #allok = TRUE

#Hierarchical clustering
htree <- hclust(dist(t(data)), method = "average")
pdf("~/desktop/HRCP.pdf", width = 100, height = 50)
plot(htree)
dev.off()

trim <- c("TCGA.49.4486.01A",
          "TCGA.78.7633.01A",
          "TCGA.73.4675.01A",
          "TCGA.69.8253.01A",
          "TCGA.78.7167.01A",
          "TCGA.73.4677.01A",
          "TCGA.49.4510.01A",
          "TCGA.78.7537.01A",
          "TCGA.97.8174.01A",
          "TCGA.78.7161.01A",
          "TCGA.MP.A4TH.01A",
          "TCGA.MP.A5C7.01A",
          "TCGA.78.7166.01A",
          "TCGA.78.7535.01A",
          "TCGA.NJ.A55R.01A",
          "TCGA.69.8255.01A",
          "TCGA.38.4626.01A",
          "TCGA.69.A59K.01A",
          "TCGA.05.4426.01A"
          )
data.subset <- data[,!(colnames(data) %in% trim)]
# plot traimed dataset
htree2 <- hclust(dist(t(data.subset)), method = "average")
pdf("~/desktop/HRCP2.pdf", width = 100, height = 50)
plot(htree2)
dev.off()

# pca
pca <- prcomp(t(data))
pca.data <- pca$x

pca.var.percent <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.data <- as.data.frame(pca.data)

library(ggplot2)
PCA<- ggplot(pca.data, aes(PC1 , PC2, PC3)) +
  geom_point() +
  geom_text(label = rownames(pca.data)) +
  labs(x = paste0("PC1:", pca.var.percent[1], "%"),
       y = paste0("PC2:", pca.var.percent[2], "%"))

pdf("~/desktop/PCA(v.0) — LUAD-TCGA.pdf", width = 100, height = 50)
plot(PCA)
dev.off()

trim.pca <- c("TCGA.55.8092.01A",
              "TCGA.05.4382.01A",
              "TCGA.86.6851.01A",
              "TCGA.55.6969.01A",
              "TCGA.44.7672.01A",
              "TCGA.69.8255.01A",
              "TCGA.99.8028.01A",
              "TCGA.44.3396.01A",
              "TCGA.55.6987.01A",
              "TCGA.78.7166.01A",
              "TCGA.73.4668.01A",
              "TCGA.NJ.A55R.01A",
              "TCGA.78.7535.01A",
              "TCGA.55.8203.01A",
              "TCGA.78.7537.01A",
              "TCGA.73.4675.01A",
              "TCGA.78.7161.01A",
              "TCGA.73.7498.01A",
              "TCGA.97.8174.01A",
              "TCGA.49.4486.01A",
              "TCGA.38.4626.01A",
              "TCGA.78.7633.01A",
              "TCGA.78.7167.01A",
              "TCGA.73.4677.01A",
              "TCGA.49.4510.01A",
              "TCGA.69.8253.01A",)


dim(data)
data.subset2 <- data[,!(colnames(data) %in% trim.pca)]
dim(data.subset2)

################################################################
final.trim <- c("TCGA.55.8092.01A",
              "TCGA.05.4382.01A",
              "TCGA.86.6851.01A",
              "TCGA.55.6969.01A",
              "TCGA.44.7672.01A",
              "TCGA.69.8255.01A",
              "TCGA.99.8028.01A",
              "TCGA.44.3396.01A",
              "TCGA.55.6987.01A",
              "TCGA.78.7166.01A",
              "TCGA.73.4668.01A",
              "TCGA.NJ.A55R.01A",
              "TCGA.78.7535.01A",
              "TCGA.55.8203.01A",
              "TCGA.78.7537.01A",
              "TCGA.73.4675.01A",
              "TCGA.78.7161.01A",
              "TCGA.73.7498.01A",
              "TCGA.97.8174.01A",
              "TCGA.49.4486.01A",
              "TCGA.38.4626.01A",
              "TCGA.78.7633.01A",
              "TCGA.78.7167.01A",
              "TCGA.73.4677.01A",
              "TCGA.49.4510.01A",
              "TCGA.69.8253.01A",
              "TCGA.49.4486.01A",
              "TCGA.78.7633.01A",
              "TCGA.73.4675.01A",
              "TCGA.69.8253.01A",
              "TCGA.78.7167.01A",
              "TCGA.73.4677.01A",
              "TCGA.49.4510.01A",
              "TCGA.78.7537.01A",
              "TCGA.97.8174.01A",
              "TCGA.78.7161.01A",
              "TCGA.MP.A4TH.01A",
              "TCGA.MP.A5C7.01A",
              "TCGA.78.7166.01A",
              "TCGA.78.7535.01A",
              "TCGA.NJ.A55R.01A",
              "TCGA.69.8255.01A",
              "TCGA.38.4626.01A",
              "TCGA.69.A59K.01A",
              "TCGA.05.4426.01A")

data.subset <- data[,!(colnames(data) %in% final.trim)]

### Normalizing ----------------------------
# Clinical data
colData <- data.frame(sample.ID = colnames(data.subset))
dds <- DESeqDataSetFromMatrix(countData =  data.subset, colData = colData, design = ~1) #not Specifying model

#suggested by WGCNA on RNASeq FAQ
dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75) #26904 genes

#Variance stabilization
dds_norm <- vst(dds75)


## WGCNA needs transposed matrix which `colnames` will be genes and `rownames` will be the samples
#Normalize counts
norm_counts <- assay(dds_norm) %>% 
  t()


#========================
#Chapter Three
#Network Construction

#Choose a set of soft thershold powers
power <- c(c(1:10) , seq(from = 12, to = 50, by = 2))

library(doParallel)
sft <- pickSoftThreshold(norm_counts, powerVector = power, networkType = "signed", verbose = 5)
# based on this matrix we will able to chose the best varablie in order to make our scale free topology
sft_Data <- sft$fitIndices

# visualization
library(ggplot2)
a1 <- ggplot(sft_Data, aes(Power, SFT.R.sq, label = Power)) +geom_point () + geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') + labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') + theme_classic()


a2 <- ggplot(sft_Data, aes(Power, mean.k., label = Power)) + geom_point () + geom_text(nudge_y = 0.1) +
  labs(x = "Power", y =  "Mean Connectivity") + theme_classic ()

pdf("~/desktop/Scale free topology model fit, signed R^2.pdf", width = 15, height = 7.5) # should use 4
a1
dev.off()

pdf("~/desktop/Mean Connectivity.pdf", width = 15, height = 7.5)
a2
dev.off()


# Convert matrix to numeric
norm_counts[] <- sapply(norm_counts, as.numeric)
soft_power = 4
temp_cor <- cor
cor <- WGCNA::cor

#memory estimate w.r.t blocksize

bwnet <- blockwiseModules(norm_counts,
                          maxBlockSize=8000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

cor <- temp_cor


# Learning section — Deprogram ___:

den.merge <- bwnet[["dendrograms"]][[1]][["merge"]]
den.height <- bwnet[["dendrograms"]][[1]][["height"]]
den.order <- bwnet[["dendrograms"]][[1]][["order"]]

write.table(den.merge,"~/desktop/den.merge (LUAD-TCGA).csv" ,  quote = F, col.names = T, row.names = T, sep = "\t" )
write.table(den.height,"~/desktop/den.height (LUAD-TCGA).csv" ,  quote = F, col.names = T, row.names = T, sep = "\t" )
write.table(den.order,"~/desktop/den.order (LUAD-TCGA).csv" ,  quote = F, col.names = T, row.names = T, sep = "\t" )

# 5. Module Eigengenes
module_eigengenes <- bwnet$MEs
# Print out a preview
head(module_eigengenes)
write.table(module_eigengenes,"~/desktop/module_eigengenes (LUAD-TCGA).csv" ,  quote = F, col.names = T, row.names = T, sep = "\t" )



# get number of genes for each module
colors <- table(bwnet$colors)
write.table(colors,"~/desktop/colors (LUAD-TCGA).csv" ,  quote = F, col.names = T, row.names = T, sep = "\t" )

# Plot the dendrogram and the module colors before and after merging underneath
pdf("~/Desktop/LUAD-TCGA-WGCNA-Dendrogram.pdf", width = 100, height = 60)
Dendogram <- plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                                 c("unmerged", "merged"), 
                                 dendroLabels = F, 
                                 addGuide = T, 
                                 hang = 0.03, 
                                 guideHang = 0.05)


#========================
#Chapter Four
#Downstream analysis

library(DEGreport)
# we are referring to `cds` which was our DEGs module in line 44
Counts <- counts(cds, normalized = T)
Counts2 <- counts(cds)
design <- as.data.frame(colData(cds))

#Size factor QC

Size.factor.QC.Plot.Normalized <- degCheckFactors(Counts[, 1:105])
Size.factor.QC.Plot <- degCheckFactors(Counts2[, 1:105])

Mean.Variance.QC.plots.normalized <- degQC(Counts, design[["group"]], pvalue = res[["pvalue"]])
Mean.Variance.QC.plots <- degQC(Counts2, design[["group"]], pvalue = res[["pvalue"]])



resCov <- degCovariates(log2(counts(cds)+0.5), 
                        colData(cds))

rownames(colData) <- colData$group
cor <- degCorCov(colData(cds))




