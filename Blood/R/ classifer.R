library(tidyverse)
library(limma)
library(oligo)
library(pheatmap)
library(caret)
library(RColorBrewer)
library(SummarizedExperiment)
library(Rtsne)
library(PCAtools)


# Add path/location functions
source(here::here("Blood/R/setup.R"))

###################### Microarray data
# All data normalized
gse <- readRDS(OUTPUT("2021-08-02_norm_GSE100153.RDS"))
ID <- "GSE100153"
###################### Meta-analysis data
# Import meta-analysis results
meta <- read_csv(OUTPUT("2021-03-01_Step6_meta-analysis.csv"))
meta$X1 <- meta[,1][[1]]

# Filter genes
qval_cutoff <- 0.05

genesqval <- meta %>% 
  #filter(meta$qval.random <= qval_cutoff ) %>% 
  filter(meta$qval.random <= qval_cutoff & I2 < 0.40 & num_exp >= 4) %>% 
  arrange(TE.random)

###################### Select genes
# Select genes in the microarray using meta-analysis results

anno_summarized <- fData(gse)  %>% 
  group_by(GENE.SYMBOL) %>% 
  dplyr::summarize(ID_first = dplyr::first(ID)) %>% 
  filter(GENE.SYMBOL %in% genesqval$X1)

feature_selection <- which(fData(gse)$ID %in% anno_summarized$ID_first)
samples_selection <- which(pData(gse)$DISEASE %in% c("CONTROL", "COPD"))

sub_gse <- gse[feature_selection,samples_selection]


###################### Prepare data for classification
df <- exprs(sub_gse)
rownames(df) <- fData(sub_gse)$GENE.SYMBOL

info <- pData(sub_gse)

anno <- data.frame(DISEASE = info$DISEASE)
rownames(anno) <- rownames(info)

# Merge in one dataframe genes and disease info
datos <- merge(t(df),anno, by = 0)
rownames(datos) <- datos$Row.names
datos <- datos[,-which(colnames(datos) == "Row.names")]

###################### Pre-heatmap  
annotation_colors = list(DISEASE = c(CONTROL = "#66CC99",COPD = "#9999CC")
                         #CLUSTER = c(Control = "#ffffff", Group1 = "#b6dbff", Group2 = "#006ddb", Group3 = "#490092", Group4 = "#009292")
                         )

paletteLength <-  30
color <- rev(colorRampPalette(c("#a30014", "white", "#004b96"))(paletteLength))
mybreaks <- myBreaks(df)

pdf(FIG(c(TODAY, "_Heatmap_", dim(df)[1], "x", dim(df)[2],"_", ID, ".pdf")),
    width = 6, height = 6)

pheatmap(df, 
         annotation = anno,
         breaks = mybreaks,
         color = color,
         cutree_cols = 5,
         show_colnames = F,
         show_rownames = F,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D",
         annotation_colors = annotation_colors,
         scale = "row",
         main = str_c("Selected genes: n = ",dim(df)[1], ", samples = ",dim(df)[2])
)
dev.off()


###################### PCA
p <- pca(df, metadata = anno)
screeplot(p)
# biplot(p,
#        colby = 'DISEASE',
#        #colkey = c(CONTROL = "#66CC99",COPD = "#9999CC"),
#        hline = 0, vline = 0,
#        title = "PCA",
#        legendPosition = 'right',
#        lab =NULL)

d_pca = as.data.frame(p$rotated[,1:3])
d_pca$DISEASE <- as.factor(p$metadata$DISEASE)
#d_pca$CLUSTER <- as.factor(annoCluster$CLUSTER)

pca <- ggplot(d_pca, aes(x=PC1, y=PC2, color = DISEASE)) +
  geom_point(size=3, alpha = 0.9) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("PC1") + ylab("PC2") +
  ggtitle("PCA") +
  theme_test(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_manual(values = c("#66CC99","#9999CC")) #c("#eeeeec","#b6dbff", "#006ddb", "#490092", "#009292")


pcas <- pairsplot(p,
          colby = 'DISEASE',
          colkey = c(CONTROL = "#66CC99",COPD = "#9999CC"))

###################### tSNE
set.seed(123)
tsne_model <- Rtsne(t(df), metadata = anno, pca=T, perplexity=10,  dims=3)
d_tsne = as.data.frame(tsne_model$Y)

d_tsne$DISEASE <- as.factor(anno$DISEASE)
rownames(d_tsne) <- rownames(t(df))

tsne <- ggplot(d_tsne, aes(x=V1, y=V2, color = DISEASE)) +
  geom_point(size=3, alpha = 0.9) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("DIM1") + ylab("DIM2") +
  ggtitle("t-SNE") +
  theme_test(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_manual(values = c("#66CC99","#9999CC")) #c("#eeeeec","#b6dbff", "#006ddb", "#490092", "#009292")

dim1 <- ggpubr::ggarrange(tsne, pca, 
                          ncol = 2, 
                          labels = c("A", "B"),
                          common.legend = TRUE)

### save plots
pdf(FIG(c(TODAY, "_DimRed_", dim(df)[1], "x", dim(df)[2],"_", ID, ".pdf")),
    width = 7, height = 10)

ggpubr::ggarrange(dim1, 
                  pcas, 
                  labels = c("","C"),
                  heights = c(1,2),
                  nrow = 2)


dev.off()

###################### Data
set.seed(234) #111
# Se crean los índices de las observaciones de entrenamiento
train <- createDataPartition(y = datos$DISEASE, p = 0.7, list = FALSE, times = 1)
TrainingSet <- datos[train, ]
TestingSet  <- datos[-train, ]

TrainingSet$DISEASE <- as.factor(TrainingSet$DISEASE)
table(TrainingSet$DISEASE)

TestingSet$DISEASE <- as.factor(TestingSet$DISEASE)
table(TestingSet$DISEASE)

#library(skimr)
#skimmed <- skim_to_wide(TrainingSet)

###################### Feature Selection
# 
# library(doMC)
# registerDoMC(cores = 4)
# 
# # Tamaño de los conjuntos de predictores analizados
# subsets <- seq(1,ncol(TrainingSet),by = 10)
# # Número de resamples para el proceso de bootstrapping
# repeticiones <- 30
# 
# ctrl_rfe <- rfeControl(functions = rfFuncs, method = "boot", number = repeticiones,
#                        returnResamp = "all", allowParallel = TRUE, verbose = FALSE)
# 
# rf_rfe <- rfe(DISEASE ~ ., data = TrainingSet,
#               sizes = subsets,
#               metric = "ROC",
#               rfeControl = ctrl_rfe,
#               ntree = 500)
# rf_rfe
# #saveRDS(rf_rfe, "Blood/output_data/rf_rfe_701g.RDS")
# #rf_rfe <- readRDS("Blood/output_data/rf_rfe_701g.RDS")
# rf_rfe$fit
# plot(rf_rfe, type = c("g", "o"))
# 
# ggplot(data = rf_rfe$results, aes(x = Variables, y = Accuracy)) +
#   geom_line() +
#   scale_x_continuous(breaks  = unique(rf_rfe$results$Variables)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = Accuracy - AccuracySD, ymax = Accuracy + AccuracySD),
#                 width = 0.2) +
#   geom_point(data = rf_rfe$results %>% slice(which.max(Accuracy)),
#              color = "red") +
#   theme_bw()

###################### Select important features
# TrainingSet <- TrainingSet[,which(colnames(TrainingSet) %in% c("DISEASE",gsub("`","",rf_rfe$optVariables)))]
# TestingSet  <- TestingSet[ ,which(colnames(TestingSet) %in% c("DISEASE",gsub("`","",rf_rfe$optVariables)))]

###################### Classifier
# SVM model (polynomial kernel)
# Build Training model
set.seed(111)

Model <- train(DISEASE ~ ., data = TrainingSet,
               method = "svmPoly",
               na.action = na.omit,
               metric='ROC',
               preProcess=c("scale","center"),
               trControl= trainControl(method="repeatedcv", 
                                       number = 10, 
                                       repeats = 3,
                                       #summaryFunction = twoClassSummary,
                                       savePredictions = T,
                                       classProbs = T),
               tuneGrid = data.frame(degree= 1,scale=0.5,C=seq(0,2, length = 10))
)

pdf(FIG(c(TODAY, "_Cost_", ID, ".pdf")), width = 6, height = 4)
plot(Model)
dev.off()
Model

# saveRDS(Model, "Blood/output_data/Model_827g.RDS")
# Model <- readRDS("Blood/output_data/Model_827g.RDS")

# Feature importance
Importance <- varImp(Model)
Importance$importance <- Importance$importance[which(Importance$importance$COPD > 70),]

pdf(FIG(c(TODAY, "_ImportanceGenes_", ID, ".pdf")), width = 4, height = 6)
plot(Importance, csi = 1)
dev.off()

Importance <- varImp(Model)

# Apply model for prediction
Model.training <-predict(Model, TrainingSet) # Apply model to make prediction on Training set
Model.testing <-predict(Model, TestingSet) # Apply model to make prediction on Testing set

# Model performance (Displays confusion matrix and statistics)
Model.training.confusion <-confusionMatrix(Model.training, TrainingSet$DISEASE)
Model.testing.confusion <-confusionMatrix(Model.testing, TestingSet$DISEASE)

print(Model.training.confusion)
print(Model.testing.confusion)

library(MLeval)
pred <- predict(Model, newdata=TestingSet, type="prob")
evalm(data.frame(pred, TestingSet$DISEASE))

svn <- rownames(Importance$importance)[which(Importance$importance$COPD > 70)]
#svn <- fData(sub_gse)$GENE.SYMBOL
#svn <- gsub("`","",rf_rfe$optVariables)

###################### Select important features
svn_df <- df[which(rownames(df) %in% svn),]

info <- pData(sub_gse)
svn_anno <- data.frame(DISEASE = info$DISEASE
                   #CLUSTER = annoCluster$CLUSTER,
                   #AGE = as.numeric(info$AGE),
                   #SEX = as.numeric(info$SEX),
                   #FEV1 = as.numeric(info$`fev1:ch1`),
                   #FEV1FVC = as.numeric(info$`fev1fvc:ch1`)
                   )
rownames(svn_anno) <- rownames(info)

# Merge in one dataframe genes and disease info
svn_datos <- merge(svn_df,svn_anno)


###################### heatmap
paletteLength <-  30
color <- rev(colorRampPalette(c("#a30014", "white", "#004b96"))(paletteLength))
mybreaks <- myBreaks(svn_df)

pdf(FIG(c(TODAY, "_Heatmap_", dim(svn_df)[1], "x", dim(svn_df)[2],"_", ID, ".pdf")), 
    width = 6, height = 6)
pheatmap(svn_df, 
         annotation = svn_anno,
         breaks=mybreaks,
         border_color = NA,
         color = color,
         cutree_cols = 5,
         show_colnames = F,
         show_rownames = T,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D",
         annotation_colors = annotation_colors,
         scale = "row",
         main = str_c("Top SVM genes: n = ",dim(svn_df)[1], ", samples = ",dim(svn_df)[2])
)
dev.off()


pdf(FIG(c(TODAY, "_BoxplotImportance_", dim(svn_df)[1],"_", ID, ".pdf")), 
    width = 15, height = 8)
featurePlot(x = TrainingSet[, svn], 
            y = TrainingSet$DISEASE, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

dev.off()

# ###################### tSNE
# set.seed(123)
# tsne_model <- Rtsne(t(svn_df), metadata = svn_anno, pca=T, perplexity=10,  dims=3)
# d_tsne = as.data.frame(tsne_model$Y)
# 
# d_tsne$DISEASE <- as.factor(annoCluster$DISEASE)
# d_tsne$CLUSTER <- as.factor(annoCluster$CLUSTER)
# rownames(d_tsne) <- rownames(t(df))
# #d_tsne_1$DISEASE <- as.factor(copdanno$DISEASE)
# 
# ggplot(d_tsne, aes(x=V1, y=V2, color = DISEASE)) +
#   geom_point(size=3, alpha = 0.9) +
#   guides(colour=guide_legend(override.aes=list(size=6))) +
#   xlab("Dim 1") + ylab("Dim 2") +
#   ggtitle("t-SNE") +
#   theme_test(base_size=20) +
#   theme(axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   scale_colour_manual(values = c("#66CC99","#9999CC")) #c("#eeeeec","#b6dbff", "#006ddb", "#490092", "#009292")
# 
# p <- pca(svn_df, metadata = svn_anno)
# screeplot(p)
# # biplot(p,
# #        colby = 'DISEASE',
# #        #colkey = c(CONTROL = "#66CC99",COPD = "#9999CC"),
# #        hline = 0, vline = 0,
# #        title = "PCA",
# #        legendPosition = 'right',
# #        lab =NULL)
# 
# d_pca = as.data.frame(p$rotated[,1:4])
# d_pca$DISEASE <- as.factor(p$metadata$DISEASE)
# d_pca$CLUSTER <- as.factor(annoCluster$CLUSTER)
# 
# ggplot(d_pca, aes(x=PC1, y=PC3, color = DISEASE)) +
#   geom_point(size=3, alpha = 0.9) +
#   guides(colour=guide_legend(override.aes=list(size=6))) +
#   xlab("PC1") + ylab("PC3") +
#   ggtitle("PCA") +
#   theme_test(base_size=20) +
#   theme(axis.text.x=element_blank(),
#         axis.text.y=element_blank()) +
#   scale_colour_manual(values = c("#66CC99","#9999CC")) #c("#eeeeec","#b6dbff", "#006ddb", "#490092", "#009292")
# 
# 
# pairsplot(p,
#           colby = 'DISEASE',
#           colkey = c(CONTROL = "#66CC99",COPD = "#9999CC"))
# 
# 

################################ Check with meta-analysis datasets
svn <- fData(sub_gse)$GENE.SYMBOL
geo <- readRDS(OUTPUT("2020-11-25_Step3_LungTissue-CURATED.RDS"))

feature_heatmap <- function(i,svn, eval_model = F){
  gse1 <- geo[[i]]
  
if(!(i %in% c(0))){


feature_selection1 <- which(fData(gse1)$GENE.SYMBOL %in% svn)
samples_selection1 <- which(pData(gse1)$DISEASE %in% c("CONTROL", "COPD"))

sub_gse1 <- gse1[feature_selection1,samples_selection1]

### 

df1 <- exprs(sub_gse1)
rownames(df1) <- fData(sub_gse1)$GENE.SYMBOL
df1 <- df1[unique(rownames(df1)),]

info1 <- as.data.frame(pData(sub_gse1))

} else {
    feature_selection1 <- which(rowData(gse1)$GENE.SYMBOL %in% svn)
    #feature_selection4 <- which(rowData(gse4)$GENE.SYMBOL %in% genesqval$X1)
    samples_selection1 <- which(colData(gse1)$DISEASE %in% c("CONTROL", "COPD"))
    
    sub_gse1 <- gse1[feature_selection1,samples_selection1]
    
    ### 
    
    
    df1 <- assay(sub_gse1)
    df1 <- log(df1+2)
    rownames(df1) <- rowData(sub_gse1)$GENE.SYMBOL
    
    info1 <- as.data.frame(colData(sub_gse1))
    
}

anno1 <- data.frame(DISEASE = info1$DISEASE)
rownames(anno1) <- rownames(info1)

paletteLength <-  30
color <- rev(colorRampPalette(c("#a30014", "white", "#004b96"))(paletteLength))
mybreaks <- myBreaks(df1)

f <- pheatmap(na.omit(df1), 
         annotation = anno1,
         breaks = mybreaks,
         annotation_colors = list(DISEASE = c(CONTROL = "#66CC99",COPD = "#9999CC")),
         color = color,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D",
         cluster_rows = T,
         show_colnames = F,
         show_rownames = F,
         cutree_cols = 3,
         scale = "row",
         main = str_c(gsub("[.].*","",names(geo))[i]," Selected genes: n = ",dim(df1)[1], ", samples = ",dim(df1)[2])
)

return(ggplotify::as.ggplot(f))

  if(eval_model == T) {
    datos1 <- merge(t(df1),anno1, by = 0)
    rownames(datos1) <- datos1$Row.names
    datos1 <- datos1[,-which(colnames(datos1) == "Row.names")]
    datos1$DISEASE <- as.factor(datos1$DISEASE)
    Model.testing <-predict(Model, datos1) 
    Model.testing.confusion <-confusionMatrix(Model.testing, datos1$DISEASE)
    print(Model.testing.confusion)
  }

}


f <- lapply(1:5, feature_heatmap, svn= svn)

pdf(FIG(c(TODAY, "_allHeatmap_", length(svn), ".pdf")), 
    width = 20, height = 15)
gridExtra::marrangeGrob(f, nrow = 2, ncol = 3)
dev.off()

######
feature_heatmap_rnaseq <- function(i, svn) {
  
  gse4 <- geo[[i]][[1]]
  
  feature_selection4 <- which(rowData(gse4)$GENE.SYMBOL %in% svn)
  #feature_selection4 <- which(rowData(gse4)$GENE.SYMBOL %in% genesqval$X1)
  samples_selection4 <- which(colData(gse4)$DISEASE %in% c("CONTROL", "COPD"))
  
  sub_gse4 <- gse4[feature_selection4,samples_selection4]
  
  ### 
  
  
  df4 <- assay(sub_gse4)
  rownames(df4) <- rowData(sub_gse4)$GENE.SYMBOL
  
  info4 <- as.data.frame(colData(sub_gse4))
  
  anno4 <- data.frame(DISEASE = info4$DISEASE)
  rownames(anno4) <- rownames(info4)
  
  # Merge in one dataframe genes and disease info
  #datos <- merge(as.data.frame(t(df)),anno)
  datos4 <- as.data.frame(t(df4))
  datos4$DISEASE <- anno4$DISEASE
  datos4 <- datos[,unique(colnames(datos4))]
  
  
  paletteLength <-  30
  color <- rev(colorRampPalette(c("#a30014", "white", "#004b96"))(paletteLength))
  mybreaks <- myBreaks(log(df4+2))
  
  f <- pheatmap(log(df4+2), 
                annotation = anno4,
                breaks = mybreaks,
                color = color,
                cluster_rows = F,
                annotation_colors = list(DISEASE = c(CONTROL = "#66CC99",COPD = "#9999CC")),
                cutree_cols = 3,
                show_colnames = F,
                show_rownames = T,
                clustering_distance_rows = "euclidean",
                clustering_method = "ward.D",
                scale = "row",
                main = str_c(names(geo)[i]," Selected genes: n = ",dim(df4)[1], ", samples = ",dim(df4)[2])
  )
  
  return(ggplotify::as.ggplot(f))
}
