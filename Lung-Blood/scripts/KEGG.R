library(tidyverse)
library(pathview)

library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)

#############
## meta analysis resulrs
blood_meta <- read_csv(here::here("Blood/output_data/Blood_2021-03-01_Step6_meta-analysis.csv")) %>% 
  filter(num_exp > 1)
lung_meta <- read_csv(here::here("Tissue/output_data/2021-03-01_Step6_meta-analysis.csv"))


# filter genes with a qvalue
qval_cutoff <- 0.05
blood <- filter(blood_meta,qval.random <= qval_cutoff & I2 < 0.40 & num_exp > 3) %>% rename(genes = X1)
lung <- filter(lung_meta,qval.random <= qval_cutoff & I2 < 0.40 & num_exp >= 9) %>% rename(genes = X1)

library(DOSE)
library(org.Hs.eg.db)
library(enrichplot)


signature <- lung
signature <- signature[,c(1,2)]
logfc <-  signature$TE.random
names(logfc) <- signature$X1

eid <- ensembldb::select(org.Hs.eg.db, names(logfc), "ENTREZID", "SYMBOL")[["ENTREZID"]]
remove <- is.na(eid)
logfc <- logfc[!remove]
names(logfc) <- eid[!remove]

logfc <- logfc[!duplicated(logfc)]

library(pathview)
library(clusterProfiler)

kk <- enrichKEGG(gene         = names(logfc),
                 organism     = 'hsa',
                 pvalueCutoff = 0.5)
head(kk)
barplot(kk)

kk2 <- gseKEGG(geneList     = sort(na.omit(logfc),decreasing = T),
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 10,
               pvalueCutoff = 0.05,
               verbose      = FALSE)
head(kk2)
dotplot(kk2) + ggtitle("dotplot for GSEA")
ridgeplot(kk2)
gseaplot2(kk2, geneSetID = 1, title = kk2$Description[1],pvalue_table = TRUE)

pathways <- c("hsa04060", "hsa04064", "hsa05323", #3
              "hsa04621", "hsa04145", "hsa05171", #6
              "hsa04514", "hsa04613", "hsa05208", #9
              "hsa04672", "hsa04061", "hsa05321", #12
              "hsa00190", "hsa05014"
              )

pathview(gene.data  = logfc*10,
         pathway.id = pathways[14],
         kegg.dir   =  here::here("Lung-Blood/output_data/kegg"),
         species    = "hsa",
         out.suffix = "lung",
         limit      = list(gene=max(abs(na.omit(logfc))), cpd=1))
