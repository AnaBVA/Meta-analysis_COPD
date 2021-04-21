library(tidyverse)
library(ggpubr)
library(oligo)

## Re shape data 
pivot_expr <- function(exprSet,GSE = 1){
  
  # identify COPD samples
  copd <- data.frame(exprs(exprSet)[,which(pData(exprSet)$DISEASE == "COPD")])
  copd <- tibble::rownames_to_column(copd, "genes")
  copd$genes <- fData(exprSet)$GENE.SYMBOL
  pes <- pivot_longer(copd,cols = -1,names_to = "GSM", values_to = "expr") %>%
    add_column(Disease = "COPD")
  
  # identify Control samples
  h <- data.frame(exprs(exprSet)[,which(pData(exprSet)$DISEASE == "CONTROL")])
  h <- tibble::rownames_to_column(h, "genes")
  h$genes <- fData(exprSet)$GENE.SYMBOL
  ph <- pivot_longer(h,cols = -1,names_to = "GSM", values_to = "expr") %>%
    add_column(Disease = "Control")
  
  # merge data as long table
  df <- rbind(pes,ph)
  df$Disease <- as.factor(df$Disease)
  df$GSE <- GSE
  df$expr <- as.numeric(df$expr)
  
  return(df)
}

##############################################################
# Lung Data
tissue <- readRDS(here::here("Tissue/output_data/2021-03-30_Step3_LungTissue-CURATED.RDS"))
#fc_tissue <- read_csv("Tissue/output_data/2020-10-07_Step5_Full_Tables.csv")


df1 <- lapply(c(1:length(tissue))[-c(4,10,12)], # 4 and 12 are RNAseq, 10 is low quality
              function(i){pivot_expr(tissue[[i]][[1]],GSE = names(tissue[i]))}
)

tdf <- bind_rows(df1)

texpplot <- ggplot(tdf,aes(x=GSE,y=expr, fill=Disease)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  theme_classic()

tdfvar <- tdf %>% 
  group_by(genes,Disease,GSE) %>%
  summarise(var=var(expr), mean=mean(expr)) 

tvarplot <- ggplot(tdfvar,aes(x=GSE,y=var, fill=Disease)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  theme_classic()

tlogvarplot <- ggplot(tdfvar,aes(x=GSE,y=log10(var), fill=Disease)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  theme_classic()

tplot <- ggarrange(texpplot,tvarplot,tlogvarplot,
                   common.legend = TRUE,
                   nrow = 3)

##############################################################
## Blood
blood <- readRDS(here::here("Blood/output_data/Blood_2020-11-25_Step3_LungTissue-CURATED.RDS"))
names(blood) <- c("GSE112811","GSE56766","GSE42057","GSE94916","GSE76705")


df2 <- lapply(c(1:length(blood)), 
              function(i){pivot_expr(blood[[i]],GSE = names(blood[i]))}
)

bdf <- bind_rows(df2)

bexpplot <- ggplot(bdf,aes(x=GSE,y=expr, fill=Disease)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  theme_classic() 

bdfvar <- bdf %>% 
  group_by(genes,Disease,GSE) %>%
  summarise(var=var(expr), mean=mean(expr)) 

bvarplot <- ggplot(bdfvar,aes(x=GSE,y=var, fill=Disease)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  theme_classic()

 ggplot(bdfvar,aes(x=GSE,y=var, fill=Disease)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  theme_classic()

blogvarplot <- ggplot(bdfvar,aes(x=GSE,y=log10(var), fill=Disease)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  theme_classic()

bplot <- ggarrange(bexpplot,bvarplot,blogvarplot,
          common.legend = TRUE,
          nrow = 3)

##############################################################
## #80624a lung
## #c45a5a blood

ggarrange(bplot,tplot,
          widths = c(1,1.5),
          ncol = 2)


## summarizing var
# bdfvar %>% select(genes,GSE, var) %>% pivot_wider(names_from = GSE, values_from = var)
bvarsum <- bdfvar %>%
  group_by(Disease,genes) %>% 
  summarise(median_var = median(var),
            mean_var = mean(var),
            tissue = "Blood")

tvarsum <- tdfvar %>%
  group_by(Disease, genes) %>% 
  summarise(median_var = median(var),
            mean_var = mean(var),
            tissue = "Lung")

varsum <- rbind(bvarsum,tvarsum)

stat.test <- varsum %>%
  group_by(Disease) %>%
  wilcox_test(median_var ~ tissue) %>%
  adjust_pvalue(method = "bonferroni") %>%
  add_significance("p.adj")

stat.test

disease_violoin <- ggplot(varsum,aes(x=Disease,y=log10(median_var), fill=tissue)) + 
  geom_violin() + 
  scale_fill_manual(values=c("#c45a75","#5abdc4")) + 
  #stat_compare_means(paired = T) +
  theme_classic(13)

tissue_violin <- ggplot(varsum,aes(x=tissue,y=log10(median_var), fill=Disease)) + 
  geom_violin() + 
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  #stat_compare_means(paired = T) +

ggarrange(disease_violin, tissue_violin,
          ncol = 2)

disease_df <- varsum %>% 
  select(genes,Disease, median_var, tissue) %>% 
  pivot_wider(names_from = Disease, values_from = median_var)

lung_scatter <- ggscatter(disease_df %>% filter(tissue == "Lung"), 
          x = "Control", y = "COPD",
          color = "#5abdc4", alpha = 0.5,
          add = "reg.line") + 
  labs(title = "Lung", x = "Control variance", y = "COPD varance") +
  stat_regline_equation(label.x = 1, label.y = 32) +
  stat_cor(label.x = 1, label.y = 34) +
  theme_bw(13)

blood_scatter <- ggscatter(disease_df %>% filter(tissue == "Blood"), 
                           x = "Control", y = "COPD",
                           color = "#c45a75", alpha = 0.5,
                           add = "reg.line") + 
  labs(title = "Blood", x = "Control variance", y = "COPD varance") +
  stat_regline_equation(label.x = 1, label.y = 32) +
  stat_cor(label.x = 1, label.y = 34) +
  theme_bw(13)

ggarrange(lung_scatter,blood_scatter,
          ncol = 2)

tissue_df <- varsum %>% 
  select(genes,Disease, median_var, tissue) %>% 
  pivot_wider(names_from = tissue, values_from = median_var)

control_scatter <- ggscatter(tissue_df %>% filter(Disease == "Control"), 
                           x = "Blood", y = "Lung",
                           color = "#66CC99", alpha = 0.5,
                           add = "reg.line") + 
  labs(title = "Control", x = "Blood variance", y = "Lung variance") +
  stat_regline_equation(label.x = 1, label.y = 32) +
  stat_cor(label.x = 1, label.y = 34) +
  theme_bw(13)

copd_scatter <- ggscatter(tissue_df %>% filter(Disease == "COPD"), 
          x = "Blood", y = "Lung",
          color = "#9999CC", alpha = 0.5,
          add = "reg.line") + 
  labs(title = "COPD", x = "Blood variance", y = "Lung variance") +
  stat_regline_equation(label.x = 1, label.y = 32) +
  stat_cor(label.x = 1, label.y = 34) +
  theme_bw(13)

ggarrange(lung_scatter,blood_scatter,control_scatter,copd_scatter,
          ncol = 2, nrow = 2)


ggplot(tissue_df %>% filter(Disease == "Control"), aes(x = Blood, y = Lung)) + geom_point()
ggplot(tissue_df %>% filter(Disease == "COPD"), aes(x = Blood, y = Lung)) + geom_point()
