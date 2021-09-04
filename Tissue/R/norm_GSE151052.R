## Normalization of GSE151052

library(GEOquery)
library(affy)
library(oligo)
library(limma)
library(tidyverse)
library(org.Hs.eg.db)
library(hugene10sttranscriptcluster.db)

# Get files
#getGEOSuppFiles("GSE151052", baseDir = here::here("Tissue/output_data/"))
geo <- getGEO("GSE151052", destdir = here::here("Tissue/output_data/"), )
geo <- geo[[1]]

# Read files
path <- "Tissue/output_data/GSE151052/GSE151052_RAW"
celfiles <- affy::list.celfiles(path, full.names=TRUE)
rawData <- oligo::read.celfiles(celfiles)

# Normalizing data
normData <- rma(rawData)
boxplot(normData)
hist(normData)

# Meta data samples
pData(normData) <- pData(geo)

# Annotation
anno_normData <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                       keys = (featureNames(normData)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

anno_normData <- subset(anno_normData, !is.na(SYMBOL))
anno <- anno_normData[which(anno_normData$PROBEID %in% featureNames(normData)),]
anno <- anno[!duplicated(anno$PROBEID),]
anno <- anno[!duplicated(anno$SYMBOL),]
rownames(anno) <- anno$PROBEID

normData <- normData[anno$PROBEID,]
fData(normData) <- anno
rownames(normData) <- fData(normData)$SYMBOL

# Add column of Gene symbol
fData(normData)$GENE.SYMBOL <- fData(normData)$SYMBOL
pData(normData)$DISEASE <- toupper(pData(normData)$"disease state:ch1")
pData(normData)$PATIENTID <- toupper(pData(normData)$"patient:ch1")
pData(normData)$AGE <- toupper(pData(normData)$"age:ch1")

saveRDS(normData,here::here("Tissue/output_data/Lung_2021-00-02_norm_GSE151052.RDS"))

## Contrast
#ct <- factor(anno$disease)
#design <- model.matrix(~0+ct)
## DE
#fit <- lmFit(y,design)
#fit2 <- eBayes(fit,trend=TRUE,robust=TRUE)

