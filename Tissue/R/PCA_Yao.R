### PCA with signature
###########################

source(here::here("Tissue/R/setup.R"))

### LIBRARIES
library(GEOquery)
library(tidyverse)
library(PCAtools)
library(cowplot)

###########################
### DATA
#meta <- read_csv(OUTPUT("2021-03-01_Step6_meta-analysis.csv"))

# 38 Biomarkers from Yao yangwei, et al.
biomarkers <- c("HBEGF","DIO2","CLCN3","SEPT4","FAT1","CTSE","CRIP1","ACADVL","CNTN3","UQCRQ","ASPN","ZNF786","RARRES2","BTC","FNDC1","DUSP1","C6orf145","NUTF2","TNN","COQ9","SCG5","BCHE","NR4A2","HS6ST3","SHE","C20orf111","REEP2","C19orf63","IRS2","FA2H","ACTL6A","NR4A3","DAO","VNN2","IGFL2","ZNF692","CAMK1D","HCAR2")

# Microarray normalized data 
geo <- readRDS(OUTPUT("2021-03-30_Step3_LungTissue-CURATED.RDS"))
names(geo)

# GSE47460
geo5 <- geo[[5]][[1]]

# Select 38 biomarkers
sig5 <- which(fData(geo5)[,"GENE_SYMBOL"] %in% biomarkers)

# Select data from ExpressionSet
geo5 <- geo5[sig5,]
metadata5 <- pData(geo5)[,c("disease state:ch1", "%emphysema (f-950):ch1","%predicted dlco:ch1","%predicted fev1 (post-bd):ch1","age:ch1","gold stage:ch1","ild subtype:ch1","Sex:ch1","smoker?:ch1")]

# Select data without ILD samples
samples <- which(metadata5$`disease state:ch1`!= "Interstitial lung disease")

# Select data from ExpressionSet
geo5_no_ILD <- geo5[,samples]
metadata5_no_ILD <- pData(geo5_no_ILD)[,c("disease state:ch1", "%emphysema (f-950):ch1","%predicted dlco:ch1","%predicted fev1 (post-bd):ch1","age:ch1","gold stage:ch1","ild subtype:ch1","Sex:ch1","smoker?:ch1")]

########################### PCA
# All samples
p <- pca(exprs(geo5), metadata = metadata5, removeVar = 0.1)

eigencorplot(p,metavars = c("disease state:ch1", "%emphysema (f-950):ch1","age:ch1","gold stage:ch1","Sex:ch1","smoker?:ch1"))
pairsplot(p,
          colby = 'disease state:ch1',
          components = getComponents(p, c(1:4)))

ild <- biplot(p,
              colby = 'disease state:ch1',
              #shape = "patient:ch1",
              hline = 0, vline = 0,
              lab = "",
              legendPosition = 'right')


# Samples without ILD
p_no_ILD <- pca(exprs(geo5_no_ILD), metadata = metadata5_no_ILD, removeVar = 0.1)
control <- biplot(p_no_ILD,
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

