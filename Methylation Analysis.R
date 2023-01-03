# Making CSV File
library(plyr)
setwd("~/desktop/Methylation Analysis/Data/GSE77954_RAW/")
all <- data.frame(data = dir(path = ".", pattern = "Grn.idat"))
Sentix_ID <- data.frame(Sentix_ID = sub("^\\w{10}_", "", all$data))
Sentix_ID <- data.frame(Sentix_ID = sub("........\\w$", "", Sentix_ID$Sentix_ID))
Sentix_ID <- data.frame(Sentix_ID = sub("_\\w+", "", Sentix_ID$Sentix_ID), Sentix_Position = sub("\\w*_", "", Sentix_ID$Sentix_ID) )

library(rvest)
library(dplyr)

link <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77954"
page <- read_html(link)

CT  <- page %>% html_nodes("tr:nth-child(25) tr td+ td") %>% html_text()
GSE <- page %>% html_nodes("tr:nth-child(25) tr a") %>% html_text()
Sample_Group <- c((rep("tumor", 34)), (rep("normal", 14 )))

GSE77954 <- data.frame(CT, GSE, stringsAsFactors = F)
GSE77954 <- data.frame(Sample_Name = sub("\\[\\w+\\]$", "", GSE77954$CT), Sample_Plate = GSE, Sample_Group = Sample_Group, stringsAsFactors = F)

PD <- data.frame(Sample_Name = GSE77954$Sample_Name,
                       Sample_Plate = GSE77954$Sample_Plate ,
                       Sample_Group = GSE77954$Sample_Group,
                       Pool_ID = NA,
                       Project = NA,
                       Sample_Well = NA,
                       Slide(Sentrix_ID) = Sentix_ID$Sentix_ID,
                       Array(Sentrix_Position) = Sentix_ID$Sentix_Position)
write.csv(csv_file, "~/desktop/Methylation Analysis/Data/GSE77954_RAW/csv_file.csv")

# Start the analysis
install.packages("ChAMP") #requiemnet
library(ChAMP)

setwd("~/desktop/Methylation Analysis/")
myImport <- champ.load("~/desktop/Methylation Analysis/Data/GSE77954_RAW/")
myLoad <- champ.filter()


#This is a useful function to demonstrate the distribution of your CpG list.
#If you get a DMP list, you may use this function to check your DMP’s distribution on chromosome.

CpG.GUI()
QC.GUI(beta=myLoad$beta,arraytype="450K") # We have diffrent results because `champ.filter()` function


myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
champ.SVD(beta=myNorm,pd=myLoad$pd)

myDMP  <- champ.DMP()
DMP.GUI()

myDMR <- champ.DMR()
DMR.GUI()

myGSEA <- champ.GSEA() 
