library(PCAtools)


geo <- readRDS("Tissue/output_data/2020-10-07_Step3_LungTissue-CURATED.RDS")
i <- geo[[5]][[1]]

#lung <- read_csv(here::here("Tissue/output_data/2021-03-01_Step6_meta-analysis.csv"))

# filter genes with a qvalue
bt # 15 genes

geneid <- rowData(i)$GENE.SYMBOL
names(geneid) <- rownames(rowData(i))
geneid <- geneid[which(geneid %in% bt$X1)]
ii <- i[names(geneid), which(colData(i)$DISEASE == "COPD" | colData(i)$DISEASE == "CONTROL" )] 
p <- pca(assay(ii), metadata = colData(ii), removeVar = 0.1)


geneid <- fData(i)$GENE.SYMBOL
names(geneid) <- rownames(fData(i))
geneid <- which(geneid %in% bt$X1)

ii <- i[geneid, which(pData(i)$DISEASE == "COPD" | pData(i)$DISEASE == "CONTROL") ] # & pData(ii)$DISEASE != "Interstitial lung disease"


p <- pca(exprs(ii), metadata = pData(ii), removeVar = 0.1)
screeplot(p)
biplot(p,
       colby = 'DISEASE',
       hline = 0, vline = 0,
       legendPosition = 'right',
       lab =NULL)

pairsplot(p,colby = 'DISEASE')


eigencorplot(p,
             metavars = c("geo_accession","DISEASE"),
             col = c('white', 'red2', 'darkred'),
             plotRsquared = TRUE,
             scale = TRUE)


library(Rtsne)

tsne_model_1 <- Rtsne(t(exprs(ii)), metadata = pData(ii), pca=TRUE, perplexity=30,  dims=3)
d_tsne_1 = as.data.frame(tsne_model_1$Y)


ggplot(d_tsne_1, aes(x=V1, y=V2, color = pData(ii)$DISEASE)) +
  geom_point(size=3, alpha = 0.8) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_brewer(palette = "Set2")

