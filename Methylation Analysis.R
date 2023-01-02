install.packages("ChAMP")
library(ChAMP)

setwd("~/desktop/Methylation Analysis/")
myLoad <- champ.load("~/desktop/Methylation Analysis/Data/GSE77954_RAW/")

rows <- c("Sample_Name", "Sample_Plate", "Sample_Group", "Pool_ID", "Project", "Sample_Well", "Sentix_ID", "Sentix_Position")
library(plyr)
paths <- dir(pattern = "\\Grn.idat$")
names(paths) <- basename(paths)
all <- ldply(paths, read.csv)


setwd("~/desktop/Methylation Analysis/Data/GSE77954_RAW/")
all <- dir(path = ".", pattern = "Grn.idat")
Sentix_ID <- grep(pattern = "94069", all, ignore.case = T, perl = FALSE, value = FALSE,
                  fixed = FALSE, useBytes = FALSE, invert = FALSE)
