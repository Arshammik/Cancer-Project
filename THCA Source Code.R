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
cn <- distinct(cn, gene_name, .keep_all = TRUE)
dim(cn)
# lose 1233 genes

rownames(cn) <- cn$gene_name
cn <- cn[,-1:-2]

#sample types
sampleID <- read.delim("../THCA/TCGA-LUAD/Data/gdc_sample_sheet.2022-06-11.tsv")
i1 <- sampleID$Data.Type == "Gene Expression Quantification"
sampleID <- sampleID[i1, , drop = F]
sampleID <- sampleID[,7:8]
a <- colnames(cn)
for(i in 1:4){a <- sub('\\.', '-', a)}
b <- data.frame(Sample.ID = a)
c <- merge(b, sampleID, by.x = "Sample.ID", by.y = "Sample.ID", sort = F)
colnames (cn) <- c$Sample.Type

# we reached to our tidy dataset

library(DESeq2)
gr <- factor(c$Sample.Type)
colSums(cn) #Through this we can demonstrates 
colData <- data.frame(group = gr, type = "paired-end")
cds <- DESeqDataSetFromMatrix(cn, colData, design = ~group)
cds <- DESeq(cds)
cnt <- log2(1+(counts(cds, normalize = T))) #getting normalized counts

write.table(cnt, "~/desktop/LUAD/TCGA-LUAD/Results/expression(log2+1)(cnt).csv",  quote = F, col.names = T, row.names = T, sep = "\t")

#DEGs
dif <- results(cds, c("group", "Solid Tissue Normal" , "Primary Tumor"))
write.table(dif, "~/desktop/Systematic Review/TCGA-LUAD/Results/dif.csv",  quote = F, col.names = T, row.names = T, sep = "\t")

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

########
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

data.trimmed <- data[,!(colnames(data) %in% final.trim)]

### Normalizing ----------------------------
# Clinical data


