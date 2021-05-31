library(tidyverse)
library(limma)
library(oligo)
library(pheatmap)
library(caret)
library(RColorBrewer)
library(SummarizedExperiment)



# Add path/location functions
source(here::here("Tissue/R/setup.R"))

###################### Microarray data
# All data normalized
geo <- readRDS(OUTPUT("2021-03-30_Step3_LungTissue-CURATED.RDS"))
names(geo)

# Select "GSE47460"
gse5 <- geo[[5]][[1]]

###################### Meta-analysis data
# Import meta-analysis results
meta <- read_csv("Tissue/output_data/2021-03-01_Step6_meta-analysis.csv")

# Filter genes
qval_cutoff <- 0.05

genesqval <- meta %>% 
  #filter(meta$qval.random <= qval_cutoff ) %>% 
  filter(meta$qval.random <= qval_cutoff & I2 < 0.40 & num_exp >= 9) %>% 
  arrange(TE.random)

###################### Select genes
# Select genes in the microarray using meta-analysis results
feature_selection <- which(fData(gse5)$GENE.SYMBOL %in% genesqval$X1)
samples_selection <- which(pData(gse5)$DISEASE %in% c("CONTROL", "COPD"))

sub_gse5 <- gse5[feature_selection,samples_selection]


###################### Subcluster COPD data
copd <- gse5[samples_selection,which(pData(gse5)$DISEASE =="COPD")]
copddf <- exprs(copd)
rownames(copddf) <- fData(copd)$GENE.SYMBOL

copdinfo <- pData(copd)
copdanno <- data.frame(DISEASE = copdinfo$DISEASE)
rownames(copdanno) <- rownames(copdinfo)

#### Dim Reduction
library(Rtsne)

set.seed(123)
tsne_model_1 <- Rtsne(t(copddf), metadata = copdanno, pca=T, perplexity=5,  dims=3)
d_tsne_1 = as.data.frame(tsne_model_1$Y)

km_clusters <- kmeans(x = d_tsne_1[, c("V1", "V2")], centers = 4, nstart = 50)
d_tsne_1$cluster <- as.factor(km_clusters$cluster)
rownames(d_tsne_1) <- rownames(t(copddf))
#d_tsne_1$DISEASE <- as.factor(copdanno$DISEASE)

ggplot(d_tsne_1, aes(x=V1, y=V2, color = cluster)) +
  geom_point(size=3, alpha = 0.9) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_test(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_manual(values = c("#b6dbff", "#006ddb", "#490092", "#009292"))


###################### Prepare data for classification
df <- exprs(sub_gse5)
rownames(df) <- fData(sub_gse5)$GENE.SYMBOL

info <- pData(sub_gse5)

anno <- data.frame(DISEASE = info$DISEASE)
rownames(anno) <- rownames(info)

### Cluster info
annoCluster <- merge(anno,d_tsne_1, by = 0, all.x = T)
rownames(annoCluster) <- rownames(anno)

annoCluster$CLUSTER <- NA

annoCluster[is.na(annoCluster$cluster),"CLUSTER"] <- 'Control'
annoCluster[which(annoCluster$cluster=="1"),"CLUSTER"] <- 'Group1'
annoCluster[which(annoCluster$cluster=="2"),"CLUSTER"] <- 'Group2'
annoCluster[which(annoCluster$cluster=="3"),"CLUSTER"] <- 'Group3'
annoCluster[which(annoCluster$cluster=="4"),"CLUSTER"] <- 'Group4'

annoCluster$DISEASE <- as.factor(annoCluster$DISEASE)
annoCluster$CLUSTER <- as.factor(annoCluster$CLUSTER)

annoCluster <- annoCluster[,c("DISEASE","CLUSTER")]

# anno <- data.frame(DISEASE = annoCluster$CLUSTER)
# rownames(anno) <- rownames(annoCluster)

# Merge in one dataframe genes and disease info
datos <- merge(t(df),anno, by = 0)
rownames(datos) <- datos$Row.names
datos <- datos[,-which(colnames(datos) == "Row.names")]


###################### Pre-heatmap  
color <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(15))
annotation_colors = list(DISEASE = c(CONTROL = "#66CC99",COPD = "#9999CC"),
                         CLUSTER = c(Control = "#ffffff", Group1 = "#b6dbff", Group2 = "#006ddb", Group3 = "#490092", Group4 = "#009292"))

pheatmap(df, 
         annotation = annoCluster,
         color = color,
         cutree_cols = 5,
         show_colnames = F,
         show_rownames = F,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D",
         annotation_colors = annotation_colors,
         scale = "row",
         main = str_c("Selected SVN model genes: n = ",dim(df)[1], ", samples = ",dim(df)[2])
)


###################### PCA
library(PCAtools)
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
d_pca$CLUSTER <- as.factor(annoCluster$CLUSTER)

ggplot(d_pca, aes(x=PC1, y=PC2, color = DISEASE)) +
  geom_point(size=3, alpha = 0.9) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("PCA") +
  theme_test(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_manual(values = c("#66CC99","#9999CC")) #c("#eeeeec","#b6dbff", "#006ddb", "#490092", "#009292")


pairsplot(p,
          colby = 'DISEASE',
          colkey = c(CONTROL = "#66CC99",COPD = "#9999CC"))

###################### tSNE
set.seed(123)
tsne_model <- Rtsne(t(df), metadata = annoCluster, pca=T, perplexity=10,  dims=3)
d_tsne = as.data.frame(tsne_model$Y)

d_tsne$DISEASE <- as.factor(annoCluster$DISEASE)
d_tsne$CLUSTER <- as.factor(annoCluster$CLUSTER)
rownames(d_tsne) <- rownames(t(df))
#d_tsne_1$DISEASE <- as.factor(copdanno$DISEASE)

ggplot(d_tsne, aes(x=V1, y=V2, color = DISEASE)) +
  geom_point(size=3, alpha = 0.9) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("") + ylab("") +
  ggtitle("t-SNE") +
  theme_test(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_manual(values = c("#66CC99","#9999CC")) #c("#eeeeec","#b6dbff", "#006ddb", "#490092", "#009292")



###################### Data
set.seed(234) #111
# Se crean los índices de las observaciones de entrenamiento
train <- createDataPartition(y = datos$DISEASE, p = 0.8, list = FALSE, times = 1)
TrainingSet <- datos[train, ]
TestingSet  <- datos[-train, ]

TrainingSet$DISEASE <- as.factor(TrainingSet$DISEASE)
TestingSet$DISEASE <- as.factor(TestingSet$DISEASE)

library(skimr)
#skimmed <- skim_to_wide(TrainingSet)

###################### Feature Selection

library(doMC)
registerDoMC(cores = 4)

# Tamaño de los conjuntos de predictores analizados
subsets <- seq(1,ncol(TrainingSet),by = 10)

# Número de resamples para el proceso de bootstrapping
repeticiones <- 30

set.seed(111)

ctrl_rfe <- rfeControl(functions = rfFuncs, method = "boot", number = repeticiones,
                       returnResamp = "all", allowParallel = TRUE, verbose = FALSE)

set.seed(234)
rf_rfe <- rfe(DISEASE ~ ., data = TrainingSet,
              sizes = subsets,
              metric = "ROC",
              rfeControl = ctrl_rfe,
              ntree = 500)

rf_rfe

#saveRDS(rf_rfe, "Tissue/output_data/rf_rfe_61g.RDS")

rf_rfe$fit

plot(rf_rfe, type = c("g", "o"))

ggplot(data = rf_rfe$results, aes(x = Variables, y = Accuracy)) +
  geom_line() +
  scale_x_continuous(breaks  = unique(rf_rfe$results$Variables)) +
  geom_point() +
  geom_errorbar(aes(ymin = Accuracy - AccuracySD, ymax = Accuracy + AccuracySD),
                width = 0.2) +
  geom_point(data = rf_rfe$results %>% slice(which.max(Accuracy)),
             color = "red") +
  theme_bw()

###################### Select important features
TrainingSet <- TrainingSet[,which(colnames(TrainingSet) %in% c("DISEASE",rf_rfe$optVariables))]
TestingSet  <- TestingSet[ ,which(colnames(TestingSet) %in% c("DISEASE",rf_rfe$optVariables))]

###################### Classifier
# SVM model (polynomial kernel)
# Build Training model
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

plot(Model)
Model

#saveRDS(Model, "Tissue/output_data/Model_61g.RDS")


# Feature importance
Importance <- varImp(Model)
plot(Importance, size = 9)

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

svn <- rownames(Importance$importance)[which(Importance$importance$COPD > 60)]
svn <- fData(sub_gse5)$GENE.SYMBOL
svn <- rf_rfe$optVariables

###################### Select important features
svn_df <- df[which(rownames(df) %in% svn),]

info <- pData(sub_gse5)
svn_anno <- data.frame(DISEASE = info$DISEASE,
                   CLUSTER = annoCluster$CLUSTER,
                   AGE = as.numeric(info$AGE),
                   EMPHYSEMA = as.numeric(info$`%emphysema (f-950):ch1`),
                   DLCO = as.numeric(info$`%predicted dlco:ch1`),
                   FEV1 = as.numeric(info$`%predicted fev1 (post-bd):ch1`),
                   GOLD = info$`gold stage:ch1`)
rownames(svn_anno) <- rownames(info)

# Merge in one dataframe genes and disease info
svn_datos <- merge(svn_df,svn_anno)


###################### heatmap
color <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(15))
pheatmap(svn_df, 
         annotation = svn_anno,
         color = color,
         cutree_cols = 5,
         show_colnames = F,
         show_rownames = T,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D",
         annotation_colors = annotation_colors,
         scale = "row",
         main = str_c("Selected SVN model genes: n = ",dim(svn_df)[1], ", samples = ",dim(svn_df)[2])
)


featurePlot(x = TrainingSet[, svn], 
            y = TrainingSet$DISEASE, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

###################### tSNE
set.seed(123)
tsne_model <- Rtsne(t(svn_df), metadata = svn_anno, pca=T, perplexity=10,  dims=3)
d_tsne = as.data.frame(tsne_model$Y)

d_tsne$DISEASE <- as.factor(annoCluster$DISEASE)
d_tsne$CLUSTER <- as.factor(annoCluster$CLUSTER)
rownames(d_tsne) <- rownames(t(df))
#d_tsne_1$DISEASE <- as.factor(copdanno$DISEASE)

ggplot(d_tsne, aes(x=V1, y=V2, color = DISEASE)) +
  geom_point(size=3, alpha = 0.9) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("Dim 1") + ylab("Dim 2") +
  ggtitle("t-SNE") +
  theme_test(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_manual(values = c("#66CC99","#9999CC")) #c("#eeeeec","#b6dbff", "#006ddb", "#490092", "#009292")

p <- pca(svn_df, metadata = svn_anno)
screeplot(p)
# biplot(p,
#        colby = 'DISEASE',
#        #colkey = c(CONTROL = "#66CC99",COPD = "#9999CC"),
#        hline = 0, vline = 0,
#        title = "PCA",
#        legendPosition = 'right',
#        lab =NULL)

d_pca = as.data.frame(p$rotated[,1:4])
d_pca$DISEASE <- as.factor(p$metadata$DISEASE)
d_pca$CLUSTER <- as.factor(annoCluster$CLUSTER)

ggplot(d_pca, aes(x=PC1, y=PC2, color = DISEASE)) +
  geom_point(size=3, alpha = 0.9) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  xlab("PC1") + ylab("PC2") +
  ggtitle("PCA") +
  theme_test(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank()) +
  scale_colour_manual(values = c("#66CC99","#9999CC")) #c("#eeeeec","#b6dbff", "#006ddb", "#490092", "#009292")


pairsplot(p,
          colby = 'DISEASE',
          colkey = c(CONTROL = "#66CC99",COPD = "#9999CC"))


###################### RandomForest
# PARALELIZACIÓN DE PROCESO
#===============================================================================
library(doMC)
registerDoMC(cores = 4)

# HIPERPARÁMETROS, NÚMERO DE REPETICIONES Y SEMILLAS PARA CADA REPETICIÓN
#===============================================================================
particiones  <- 10
repeticiones <- 5

# Hiperparámetros
hiperparametros <- expand.grid(mtry = c(3, 4, 5, 7, 9),
                               min.node.size = c(2, 3, 4, 5, 10, 15, 20, 30),
                               splitrule = "gini")

set.seed(345)
seeds <- vector(mode = "list", length = (particiones * repeticiones) + 1)
for (i in 1:(particiones * repeticiones)) {
  seeds[[i]] <- sample.int(1000, nrow(hiperparametros))
}
seeds[[(particiones * repeticiones) + 1]] <- sample.int(1000, 1)

# DEFINICIÓN DEL ENTRENAMIENTO
#===============================================================================
control_train <- trainControl(method = "repeatedcv", number = particiones,
                              repeats = repeticiones, seeds = seeds,
                              returnResamp = "final", verboseIter = FALSE,
                              allowParallel = TRUE)

# AJUSTE DEL MODELO
# ==============================================================================
set.seed(123)
modelo_rf <- train(DISEASE ~ ., data = TrainingSet,
                   method = "ranger",
                   tuneGrid = hiperparametros,
                   metric = "Accuracy",
                   trControl = control_train,
                   # Número de árboles ajustados
                   num.trees = 500)
modelo_rf
modelo_rf$finalModel

# REPRESENTACIÓN GRÁFICA
# ==============================================================================
ggplot(modelo_rf, highlight = TRUE) +
  scale_x_continuous(breaks = 1:30) +
  labs(title = "Evolución del accuracy del modelo Random Forest") +
  guides(color = guide_legend(title = "mtry"),
         shape = guide_legend(title = "mtry")) +
  theme_bw()

# Apply model for prediction
modelo_rf.training <-predict(modelo_rf, TrainingSet) # Apply model to make prediction on Training set
modelo_rf.testing <-predict(modelo_rf, TestingSet) # Apply model to make prediction on Testing set


# Model performance (Displays confusion matrix and statistics)
modelo_rf.training.confusion <-confusionMatrix(modelo_rf.training, TrainingSet$DISEASE)
modelo_rf.testing.confusion <-confusionMatrix(modelo_rf.testing, TestingSet$DISEASE)

print(modelo_rf.training.confusion)
print(modelo_rf.testing.confusion)

library(plotROC)
ggplot(modelo_rf$pred[])

################################ External validation


feature_heatmap <- function(i,svn, eval_model = F){
  gse1 <- geo[[i]]
  
  if(!(i %in% c(4,12))){
    
    
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
  
  color <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(15))
  f <- pheatmap(na.omit(df1), 
                annotation = anno1,
                annotation_colors = list(DISEASE = c(CONTROL = "#66CC99",COPD = "#9999CC")),
                color = color,
                clustering_distance_rows = "euclidean",
                clustering_method = "ward.D",
                cluster_rows = T,
                show_colnames = F,
                show_rownames = T,
                cutree_cols = 3,
                scale = "row",
                main = str_c(names(geo)[i]," Selected genes: n = ",dim(df1)[1], ", samples = ",dim(df1)[2])
  )
  
  if(eval_model == T) {
    datos1 <- merge(t(df1),anno1, by = 0)
    rownames(datos1) <- datos1$Row.names
    datos1 <- datos1[,-which(colnames(datos1) == "Row.names")]
    datos1$DISEASE <- as.factor(datos1$DISEASE)
    Model.testing <-predict(Model, datos1) 
    Model.testing.confusion <-confusionMatrix(Model.testing, datos1$DISEASE)
    print(Model.testing.confusion)
  }
  
  return(ggplotify::as.ggplot(f))
  
}


f <- sapply(1:12, feature_heatmap)
g <- gridExtra::marrangeGrob(f, nrow = 4, ncol = 3)

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
  
  
  color <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(15))
  f <- pheatmap(log(df4+2), 
                annotation = anno4,
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
