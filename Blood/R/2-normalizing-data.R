source(here::here("R/setup.R"))

library(limma)
library(oligo)
library(tidyverse)

rawCEL_normCEL <- function(gse){
  
  # select and untar files
  gsepath <- file.path(here::here("GEOquery"),gse)
  print(gse)
  
  if(length(list.files(gsepath)) == 0){
    message("No Supplementary files found")
    } else if (length(list.files(gsepath)) < 3 & length(list.files(gsepath)) != 0 ){
    ifelse(any(list.files(gsepath,pattern="tar") !=0),
                untar(list.files(gsepath,full.names=T,pattern="tar"), exdir=gsepath),
                list.files(gsepath))
    } else{
    print(length(list.files(gsepath)))
    }
  
  # for CEL files "AFFYMETRIX"
  if(any(grepl("cel",list.files(gsepath),ignore.case = T))){
    message("This is an Affymetrix array")
    celfiles <- list.celfiles(gsepath,
                              full.names=TRUE,
                              listGzipped=TRUE)
    
    # read CEL files in R
    rawData <- read.celfiles(celfiles)
    #### Figures of raw data
    #pdf(str_c("raw_",gse,"_boxplot",TODAY,".pdf"))
    ## plots of raw data
    #boxplot(rawData,target="core")
    #hist(rawData,target="core")
    #dev.off()
    ## RMA normalization
    normData <- rma(rawData)
    #### Figures of Normalized data
    #pdf(str_c("norm_",gse,"_boxplot",TODAY,".pdf"))
    ## `plots of norm data
    boxplot(normData)
    hist(normData)
    # dev.off()
    #write.csv(exprs(normData),str_c("/normalized",gse,"_normData",TODAY,".tx  t"),quote=F)
    return(normData)
  }
  
  #for TXT files "AGILENT"
  else if(sum(grepl("txt",list.files(gsepath),ignore.case = T)) >2){
    message("This is an Agilent array")
    files <- list.files(gsepath,full.names=T, pattern = "txt")
    files <- grep("annot",files,value=T,invert=T)
   
     # read agilent files
    RG <- read.maimages(files,source = "agilent.median",green.only=T)
    #plots: raw data
    boxplot(RG$E)
    hist(RG$E)
    # normalizing data
    RG <- limma::backgroundCorrect(RG, method="normexp", offset=1)
    RG$E <- normalizeBetweenArrays(RG$E, method="quantile")
    RG$E <- log2(RG$E)
    #plots: norm
    boxplot(RG$E)
    hist(RG$E,100)
    
    return(RG)
  }
  
  # anyother case
  else{
    message("This is not Affymetrix or Agilent")
  }
}



gse_table <- read_tsv(DATA("2-Table_GSE-info.txt"))

geo_norm <- sapply(gse_table$Accession,rawCEL_normCEL)
saveRDS(geo_norm,OUTPUT(c(TODAY,"_normData.xz")),compress = "xz")

sessionInfo()
