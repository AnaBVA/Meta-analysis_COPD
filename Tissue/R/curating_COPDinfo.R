library(GEOquery)
library(mgsub)
library(annotate)
library(hugene10sttranscriptcluster.db)
library(EnsDb.Hsapiens.v79)
library(tidyverse)
library(gsubfn)

### TISSUE
exp_tissue <- c('GSE27597','GSE106986','GSE37768','GSE57148','GSE47460',
                'GSE38974','GSE8581','GSE38934','GSE69818','GSE63073',
                'GSE17770','GSE1650','GSE1122','GSE119040','GSE103174','GSE124180')

######### READ data
## Experiments had been downloaded and saved in a RDS object to have an easy and
## quicker access. 
#geo <- sapply(exp_tissue, getGEO)
#saveRDS(geo,here::here("data/2020-09-03-GSE_LungTissue.RDS"))
geo <- readRDS(here::here("data/2020-09-03-GSE_LungTissue.RDS"))

######## NEW COLUMN  with homogenized vocabulary
## Here we add a new column and rename the description that each authors had
## used to describe COPD or Control samples

######## SUMMARY TABLE
## Summary information about the experiment
df <- tibble("GSE"= NA,
             "CONTROL" = NA,
             "COPD" = NA,
             "COUNTRY" = NA,
             "SUBMISSION_DATE" = NA,
             "PLATFORM" = NA,
             "NUM_GENES" = NA
                 )


########################################### 
## 1
# Select the experiment
i <- 1
p <- pData(geo[[i]][[2]])
head(p)
dis <- "copd status:ch1" 

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("0","1")

new_labels <- c("CONTROL","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)


# other infromation
p$PATIENT <- p$`patient:ch1`
p$AGE <- p$`age:ch1`
p$SEX <- p$`Sex:ch1`
p$PACKYEAR <- p$`pack years:ch1`

p$SMOKING_STATUS <- NA
p$SMOKING_STATUS[grep("non-smoker",p$title)] <- "NON-SMOKER"
p$SMOKING_STATUS[grep(" smoker",p$title)] <- "SMOKER"     

# add it to the expression object
pData(geo[[i]][[2]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[2]])
head(f)

# Parse column to have Gene.Symbol
ff <- dplyr::mutate(f,
                   geneInfo = strsplit(gene_assignment, " /// "),
                   gene1 = sapply(geneInfo, `[`, 1),
                   GENE.SYMBOL = sapply(strsplit(gene1, " // "), `[`, 2)
)

rownames(ff) <- rownames(f)
f <- dplyr::select(ff,ID,GB_LIST,GENE.SYMBOL)

length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[2]]) <- f

########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)



########################################### 
i <- 2
p <- pData(geo[[i]][[1]])
head(p)
dis <- "source_name_ch1"

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("tumor-free tissue of non-smoker",
                "tumor-free tissue of ex-smoker with COPD",
                "tumor-free tissue of smoker with COPD")

new_labels <- c("CONTROL","COPD","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)

# other infromation
p$AGE <- p$`age:ch1`
p$SEX <- p$`gender:ch1`
p$SMOKING_STATUS <- mgsub(p[,dis],old_labels,c("NON-SMOKER","EX-SMOKER","SMOKER"))

# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$GENE_SYMBOL
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f

########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)


########################################### 
## 3
# Select the experiment
i <- 3
p <- pData(geo[[i]][[1]])
head(p)
dis <- "characteristics_ch1.1"

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("phenotype: healthy smoker","phenotype: Nonsmoker","phenotype: moderate COPD")

new_labels <- c("CONTROL","CONTROL","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)

# other infromation
p$SMOKING_STATUS <- mgsub(p[,dis],old_labels,c("SMOKER","NON-SMOKER",NA))


# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$`Gene Symbol`
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f

########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)

########################################### 
library(recount)
## 4
# Select the experiment
i <- 4

## Download data from Recount2
#url <- download_study('SRP041538',outdir = file.path('download/SRP041538'))
load(file.path('download/SRP041538', 'rse_gene.Rdata'))

## Scale counts by taking into account the total coverage per sample
rse <- scale_counts(rse_gene)
rse

# saving the geoquery object 
geo[[i]][[2]] <- geo[[i]][[1]]

## adding recount2 data to the list 
geo[[i]][[1]] <- rse

### 
p <- colData(geo[[i]][[1]])
head(p)
dis <- "characteristics"

######### DISEASE information
table(unlist(p[,dis]))
ss <- sapply(p[,dis],unlist)
ss <- ss[1,]

# identify the way authors described COPD or CONTROLS
old_labels <- c("disease state: Normal","disease state: COPD")
new_labels <- c("CONTROL","COPD")

# rename 
p$DISEASE <- mgsub(ss,old_labels,new_labels)
table(p$DISEASE)

# other infromation

# add it to the expression object
colData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- rowData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- as.character(as.data.frame(f)$symbol)
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
rowData(geo[[i]][[1]]) <- f

########### Summary
pp <- pData(geo[[i]][[2]])
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(pp$contact_country),
        "SUBMISSION_DATE" = paste(unique(pp$submission_date),collapse = ","),
        "PLATFORM" = unique(pp$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)

########################################### 
## 5
# Select the experiment
i <- 5
p <- pData(geo[[i]][[1]])
head(p)
dis <- "disease state:ch1" 

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("Control",
                "Chronic Obstructive Lung Disease")
new_labels <- c("CONTROL","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)

# other infromation
p$AGE <- p$`age:ch1`
p$SMOKING_STATUS <- p$`smoker?:ch1`

# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$GENE_SYMBOL
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f


########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)


########################################### 
## 6
# Select the experiment
i <- 6
p <- pData(geo[[i]][[1]])
head(p)
dis <- "diagnosis:ch1" 

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("Normal",
                "COPD")
new_labels <- c("CONTROL","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)

# other infromation
p$AGE <- p$`age:ch1`
p$SEX <- p$`gender (1=male, 2=female):ch1`
p$PACKYEAR <- p$`pkyrs:ch1`
p$SMOKING_STATUS <- rep("SMOKER",nrow(p))

# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$GENE_SYMBOL
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f


########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)






########################################### 
## 7
# Select the experiment
i <- 7
p <- pData(geo[[i]][[1]])
head(p)

table(gsub('[0-9]+', '', p$title))
p$DISEASE <- "CONTROL"
p$DISEASE[grep("COPD",p$title)] <- "COPD"
p$DISEASE[grep("Unclassified",p$title)] <- "Unclassified"
table(p$DISEASE)


# other infromation
pp <- p$characteristics_ch1

p$AGE <- strapplyc(pp, "Age:.[0-9]*", simplify = TRUE)
p$AGE <- gsub("[^0-9.-]", "", p$AGE)

p$SEX <- strapplyc(pp, "Gender:.[A-Z]*", simplify = TRUE)
p$SEX <-  gsub("[^M|F.-]", "", p$SEX)

# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$`Gene Symbol`
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f


########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)


########################################### 
## 8
# Select the experiment
i <- 8
p <- pData(geo[[i]][[1]])
head(p)
p$MOD_title <- stringr::str_count(as.character(p$title,"L"))
dis <- "MOD_title" 

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("3",
                "2")
new_labels <- c("CONTROL","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)

# other infromation

# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$`Gene Symbol`
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f

########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)

########################################### 
## 9
# Select the experiment
i <- 9
p <- pData(geo[[i]][[1]])
head(p)
dis <- "description" 

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("Normal lung",
                "'usual' emphysemic lung",
                "Alpha-1 Antitrypsin Deficiency-related emphysemic lung")
new_labels <- c("CONTROL","COPD","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)

# other infromation

# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$`Gene Symbol`
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f


########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)

########################################### 
## 10
# Select the experiment
i <- 10
p <- pData(geo[[i]][[1]])
head(p)
dis <- "disease state:ch1" 

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("non-emphysema","emphysema")
new_labels <- c("CONTROL","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)

# other infromation

# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$`Gene Symbol`
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f


########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)


########################################### 
## 11
# Select the experiment
i <- 11
p <- pData(geo[[i]][[1]])
head(p)
dis <- "copd:ch1" 

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("no",
                "yes")
new_labels <- c("CONTROL","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)

# other infromation
p$AGE <- p$`age:ch1`
p$SEX <- p$`Sex:ch1`
p$PACKYEAR <- p$`pack years:ch1`
p$SMOKING_STATUS <- p$`smoking:ch1`

# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$`Gene Symbol`
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f


########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)

########################################### 
## 12
# Select the experiment
i <- 12
p <- pData(geo[[i]][[1]])

############ Get gene counts
#pp <- getGEOSuppFiles(names(geo)[i])
pp <- read.csv(here::here("output_data/GSE124180_gene_count_table.tsv"),sep = "\t")
rownames(pp) <- pp$ENSEMBL_GENEID
pp <- pp[,which(colnames(pp)!="X" & colnames(pp)!= "ENSEMBL_GENEID")]
colnames(pp[,p$title])

# check order of colnames 
identical(colnames(pp),p$title)

#rename colnames of gene counts df
colnames(pp) <- p$geo_accession

#Create the annotation
phenoData <- new("AnnotatedDataFrame",
                 data=pData(geo[[i]][[1]]))

# Get gene annotation
gs <- select(EnsDb.Hsapiens.v79, 
             rownames(pp),
             c("GENEID","GENENAME","SEQNAME","GENESEQSTART","GENESEQEND"), "GENEID")

rownames(gs) <- gs$GENEID
fData <- new("AnnotatedDataFrame",
             data=gs)

# Select genes with annotation
pp <- pp[which(rownames(pp) %in% gs$GENEID),]

#Create the object expressionSet
geo[[i]][[1]] <- ExpressionSet(assayData=as.matrix(pp),
                               phenoData=phenoData,
                               featureData = fData)

#### Subset data to  have only lung tissue
p <- pData(geo[[i]][[1]])
table(p$`cell type:ch1`)
geo[[i]][[1]] <- geo[[i]][[1]][,which(p$`cell type:ch1` == "bronchial epithelium")]

p <- pData(geo[[i]][[1]])
head(p)
dis <- "copd:ch1" 

######### DISEASE information
table(p[,dis])

# identify the way authors described COPD or CONTROLS
old_labels <- c("cont",
                "case")
new_labels <- c("CONTROL","COPD")

# rename 
p$DISEASE <- mgsub(p[,dis],old_labels,new_labels)
table(p$DISEASE)

# other infromation
p$AGE <- p$`age:ch1`
p$SEX <- p$`Sex:ch1`
p$PACKYEAR <- p$`pack years:ch1`
p$SMOKING_STATUS <- p$`smoke:ch1`

# add it to the expression object
pData(geo[[i]][[1]]) <- p

########### Gene.Symbol
# We also add a new column with Gene.Symbol annotation
f <- fData(geo[[i]][[1]])
head(f)
f$GENE.SYMBOL <- f$GENENAME
length(f$GENE.SYMBOL)
length(unique(f$GENE.SYMBOL))
fData(geo[[i]][[1]]) <- f


########### Summary
df <- df %>% add_row(
        "GSE" = names(geo)[i],
        "CONTROL" = table(p$DISEASE)[1],
        "COPD" = table(p$DISEASE)[2],
        "COUNTRY" = unique(p$contact_country),
        "SUBMISSION_DATE" = paste(unique(p$submission_date),collapse = ","),
        "PLATFORM" = unique(p$platform_id),
        "NUM_GENES" = length(unique(f$GENE.SYMBOL))
)




##############################################################
saveRDS(geo,here::here("data/2020-09-GSE_LungTissue-CURATED.RDS"))

df <- df[which(df$GSE != "NA"),]
write_csv(df,"output_data/2020-09-03-GSE_Summary.csv")

