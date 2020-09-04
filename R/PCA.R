### PCA with signature
###########################3

PATH = here::here()
DATA_DIR = file.path(PATH,"data")
OUTPUT_DIR = file.path(PATH,"output_data")
FIG_DIR = file.path(PATH,"fig")
TODAY = Sys.Date()

### LIBRARIES
library(GEOquery)
library(tidyverse)

### DATA
gse_table <- read.csv(file.path(DATA_DIR,"GSE_table.csv"), row.names = 1)
#geo <- sapply(rownames(gse_table), getGEO)
saveRDS(geo,here::here("data/geo_7GSEs.RDS"))


meta <- read.csv("output_data/meta-analysis_preprocessData-2020-08-13.csv")[,-1]
signature <- meta[meta$qval <= 0.01 & meta$num_exp == 7,"QEp"]
names(signature) <- meta[meta$qval <= 0.01 & meta$num_exp == 7,"genes"]

biomarkers <- c("HBEGF","DIO2","CLCN3","SEPT4","FAT1","CTSE","CRIP1","ACADVL","CNTN3","UQCRQ","ASPN","ZNF786","RARRES2","BTC","FNDC1","DUSP1","C6orf145","NUTF2","TNN","COQ9","SCG5","BCHE","NR4A2","HS6ST3","SHE","C20orf111","REEP2","C19orf63","IRS2","FA2H","ACTL6A","NR4A3","DAO","VNN2","IGFL2","ZNF692","CAMK1D","HCAR2")

head(signature)
length(signature)


library(PCAtools)

# metadata
## GSE1122
geo1 <- geo[[1]][[1]]
sig1 <- which(fData(geo1)[,"Gene Symbol"] %in% names(signature))
geo1 <- geo1[sig1,]

metadata <- pData(geo1)[,c("geo_accession","description")]
geo1 <- geo1[,which(metadata$description != "Alpha-1 Antitrypsin Deficiency-related emphysemic lung")]
metadata <- pData(geo1)[,c("geo_accession","description")]

p <- pca(exprs(geo1), metadata = metadata, removeVar = 0.1)
screeplot(p)
biplot(p,
       colby = 'description',
       hline = 0, vline = 0,
       legendPosition = 'right')

pairsplot(p,colby = 'description')
eigencorplot(p,
             metavars = c("geo_accession","description"),
             col = c('white', 'red2', 'darkred'),
             plotRsquared = TRUE,
             scale = TRUE)

## GSE1650
geo2 <- geo[[2]][[1]]
sig2 <- which(fData(geo2)[,"Gene Symbol"] %in% names(signature))
geo2 <- geo2[sig2,]

metadata2 <- pData(geo2)[,c("title","description")]
metadata2$disease <- as.factor(str_count(as.character(metadata2$title),"L"))
p <- pca(exprs(geo2), metadata = metadata2, removeVar = 0.1)
screeplot(p)
biplot(p,
       colby = 'disease',
       hline = 0, vline = 0,
       legendPosition = 'right')

eigencorplot(p,metavars = c("title",'disease'),
             col = c('white', 'red2', 'darkred'),
             plotRsquared = TRUE,
             scale = TRUE)

pairsplot(p,
          # components = getComponents(p, c(7:10)),
          colby = 'disease')


## GSE27597
geo3 <- geo[[3]][[2]]
genesymbol <- fData(geo3)[,"gene_assignment"]
gs <- vapply(strsplit(genesymbol,"//"), `[`, 2, FUN.VALUE=character(1))
gs <- gsub(" ","",gs)

sig3 <- which(gs %in% names(signature))
geo3 <- geo3[sig3,]

metadata3 <- pData(geo3)[,c("age:ch1", "copd status:ch1","lm:ch1","pack years:ch1", "patient:ch1","Sex:ch1" ,"slice:ch1")]
p <- pca(exprs(geo3), metadata = metadata3, removeVar = 0.1)
screeplot(p)
biplot(p,
       colby = 'copd status:ch1',
       #shape = "patient:ch1",
       hline = 0, vline = 0,
       legendPosition = 'right')

eigencorplot(p,col = c('white', 'red2', 'darkred'),
             plotRsquared = TRUE,
             scale = TRUE,
             metavars = c(c("age:ch1", "copd status:ch1","lm:ch1","pack years:ch1", "patient:ch1","Sex:ch1" ,"slice:ch1")))

pairsplot(p,
          #components = getComponents(p, c(1:5)),
          colby = 'copd status:ch1')

## GSE37768
geo4 <- geo[[4]][[1]]

sig4 <- which(unique(fData(geo4)[,"Gene Symbol"]) %in% names(signature))
geo4 <- geo4[sig4,]

metadata4 <- pData(geo4)[,c("phenotype:ch1", "source_name_ch1")]
p <- pca(exprs(geo4), metadata = metadata4, removeVar = 0.1)
screeplot(p)
biplot(p,
       colby = 'phenotype:ch1',
       #shape = "patient:ch1",
       hline = 0, vline = 0,
       legendPosition = 'right')

eigencorplot(p,metavars = c("phenotype:ch1", "source_name_ch1"))

pairsplot(p,
          colby = 'phenotype:ch1',
          components = getComponents(p, c(1:6)))


## GSE47460
geo5 <- geo[[5]][[1]]

sig5 <- which(fData(geo5)[,"GENE_SYMBOL"] %in% names(signature))
# sig5 <- which(fData(geo5)[,"GENE_SYMBOL"] %in% biomarkers)
geo5 <- geo5[sig5,]

metadata5 <- pData(geo5)[,c("disease state:ch1", "%emphysema (f-950):ch1","%predicted dlco:ch1","%predicted fev1 (post-bd):ch1","age:ch1","gold stage:ch1","ild subtype:ch1","Sex:ch1","smoker?:ch1")]

geo5 <- geo5[,which(metadata5$`disease state:ch1`!= "Interstitial lung disease")]
metadata5 <- pData(geo5)[,c("disease state:ch1", "%emphysema (f-950):ch1","%predicted dlco:ch1","%predicted fev1 (post-bd):ch1","age:ch1","gold stage:ch1","ild subtype:ch1","Sex:ch1","smoker?:ch1")]


p <- pca(exprs(geo5), metadata = metadata5, removeVar = 0.1)
screeplot(p)
biplot(p,
       colby = 'disease state:ch1',
       #shape = "patient:ch1",
       hline = 0, vline = 0,
       lab = "",
       legendPosition = 'right')

eigencorplot(p,metavars = c("disease state:ch1", "%emphysema (f-950):ch1","age:ch1","gold stage:ch1","Sex:ch1","smoker?:ch1"))

pairsplot(p,
          colby = 'disease state:ch1',
          components = getComponents(p, c(1:4)))

##############################
### Biomarkers 

ild <- biplot(p,
              colby = 'disease state:ch1',
              #shape = "patient:ch1",
              hline = 0, vline = 0,
              lab = "",
              legendPosition = 'right')

control <- biplot(p,
                  colby = 'disease state:ch1',
                  #shape = "patient:ch1",
                  hline = 0, vline = 0,
                  lab = "",
                  legendPosition = 'right')

pca_plots <- plot_grid(ild, control,
                     ncol = 1,
                     nrow = 2,
                     labels = c('A)', 'B)'),
                     label_fontfamily = 'serif',
                     label_fontface = 'bold',
                     label_size = 22,
                     rel_heights = c(2, 2),
                     align = 'h')

# now add the title
title <- ggdraw() + 
        draw_label(
                "38 Biomarkers publicated by Yao et.al.",
                fontface = 'bold',
                size = 22
        ) + theme(plot.margin = margin(0, 0, 0, 7))

plot_grid(
        title, pca_plots,
        ncol = 1,
        # rel_heights values control vertical title margins
        rel_heights = c(0.1, 1)
)

##############################
## GSE8581
geo7 <- geo[[7]][[1]]

sig7 <- which(fData(geo7)[,"Gene Symbol"] %in% names(signature))
geo7 <- geo7[sig7,]

metadata7 <- pData(geo7)[,c("title","characteristics_ch1")]
#f <- str_split(metadata7$characteristics_ch1,",")

metadata7 <- pData(geo7)[,c("title","characteristics_ch1")]
metadata7$disease <- factor(str_remove_all(pData(geo7)[,"title"], "[0-9]"))


p <- pca(exprs(geo7), metadata = metadata7, removeVar = 0.1)
screeplot(p)
biplot(p,
       colby = 'disease',
       #shape = "patient:ch1",
       hline = 0, vline = 0,
       legendPosition = 'right')

eigencorplot(p,metavars = c("disease","title","characteristics_ch1"))

pairsplot(p,
          colby = 'disease',
          components = getComponents(p, c(1:4)))

