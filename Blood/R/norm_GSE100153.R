## Normalization of GSE100153

library(GEOquery)
library(limma)
library(tidyverse)

# Get files
#getGEOSuppFiles("GSE100153", baseDir = here::here("Blood/output_data/"))
geo <- getGEO("GSE100153", destdir = here::here("Blood/output_data/"))
geo <- geo[[1]]

# Read files
file <- "Blood/output_data/GSE100153/GSE100153_COPD_non-normalized.txt"
x <- read.ilmn(file,expr="44876533",probeid="ID_REF")
boxplot(log2(x$E),range=0,ylab="log2 intensity")

# Meta data order
anno <- data.frame(
  id = str_extract(string = pData(geo)$title, pattern = "[0-9][0-9]_."),
  disease = pData(geo)$"sample class:ch1",
  geo = pData(geo)$"geo_accession"
)

# Normalization
y <- neqc(x)

# Merge data with ExpressionSet Object
y <- y[, anno$id]
plotMDS(y,labels=anno$disease, pch = 1 )

y2 <- as.matrix(y$E)
colnames(y2) <- anno$geo
exprs(geo) <- y2

# Filtered unexpressed probes
expressed <- rowSums(y$other$Detection < 0.05) >= 3
geo <- geo[expressed,]

boxplot(exprs(geo))

# Add column of Gene symbol
fData(geo)$GENE.SYMBOL <- fData(geo)$Symbol
pData(geo)$DISEASE <- toupper(pData(geo)$"sample class:ch1")

saveRDS(geo,here::here("Blood/output_data/Blood_2021-08-02_norm_GSE100153.RDS"))

## Contrast
#ct <- factor(anno$disease)
#design <- model.matrix(~0+ct)
## DE
#fit <- lmFit(y,design)
#fit2 <- eBayes(fit,trend=TRUE,robust=TRUE)

