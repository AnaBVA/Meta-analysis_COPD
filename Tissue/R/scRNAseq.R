library(tidyverse)
library(Seurat)

#  qrsh -l h_vmem=30G

sce <- readRDS("data/GSE173896_COPD.rds")

sce <- RunUMAP(sce, dims = 1:50)

sce$DIAGNOSIS <- "CONTROL"
sce$DIAGNOSIS[which(sce$FEV1.FVC < 70)] <- "COPD"

pdf("fig/umap.pdf")
DimPlot(sce, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()

pdf("fig/umapFEV1FVC.pdf", width = 12)
p1 <- DimPlot(sce, reduction = "umap", group.by = "DIAGNOSIS") 
p2 <- DimPlot(sce, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend()
p1 + p2
dev.off()

top = c("SELE", "SFRP2", "CXCL13", "TMEM100", "GPIHBP1", "LRRC46")

  
  
pdf("fig/umapTop.pdf", width = 10, height = 12)
FeaturePlot(sce, features = top)
dev.off()

pdf("fig/dotTop.pdf", width = 10, height = 12)
DotPlot(sce, features = top, split.by = "DIAGNOSIS") + RotatedAxis()
dev.off()

VlnPlot(pbmc, features = c("SELE", "SFRP2", "CXCL13", "TMEM100", "GPIHBP1", "LRRC46"))


# Add path/location functions
source(here::here("Tissue/R/setup.R"))

###################### Data
sce <- readRDS(DATA("GSE173896_COPD.rds"))

###################### Meta-analysis data
# Import meta-analysis results
meta <- read_csv(OUTPUT("2021-03-01_Step6_meta-analysis.csv"))
meta$X1 <- meta[,1][[1]]

# Filter genes
qval_cutoff <- 0.05

genesqval <- meta %>% 
  #filter(meta$qval.random <= qval_cutoff ) %>% 
  filter(meta$qval.random <= qval_cutoff & I2 < 0.40 & num_exp >= 9) %>% 
  arrange(TE.random)

write_csv(genesqval, OUTPUT("2021-03-01_Step6_161meta-analysis.csv"))

