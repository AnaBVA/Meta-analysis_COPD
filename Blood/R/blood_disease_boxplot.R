library(tidyverse)
library(oligo)

# Data
blood <- readRDS(here::here("Blood/output_data/Blood_2020-11-25_Step3_LungTissue-CURATED.RDS"))
#fc_blood <- read_csv("Blood/output_data/Blood_2020-11-28_Step5_Full_Tables.csv")

pivot_expr <- function(exprSet,GSE = 1){
  
  # identify COPD samples
  copd <- data.frame(exprs(exprSet)[,which(pData(exprSet)$DISEASE == "COPD")])
  copd <- tibble::rownames_to_column(copd, "genes")
  pes <- pivot_longer(copd,cols = -1,names_to = "GSM", values_to = "expr") %>%
    add_column(Disease = "COPD")
  
  # identify Control samples
  h <- data.frame(exprs(exprSet)[,which(pData(exprSet)$DISEASE == "CONTROL")])
  h <- tibble::rownames_to_column(h, "genes")
  ph <- pivot_longer(h,cols = -1,names_to = "GSM", values_to = "expr") %>%
    add_column(Disease = "Control")
  
  # merge data as long table
  df <- rbind(pes,ph)
  df$Disease <- as.factor(df$Disease)
  df$GSE <- GSE
  
  return(df)
}

df1 <- pivot_expr(blood[[1]],GSE ="GSE112811")
df2 <- pivot_expr(blood[[2]],GSE ="GSE56766")
df3 <- pivot_expr(blood[[3]],GSE ="GSE42057")
df4 <- pivot_expr(blood[[4]],GSE ="GSE94916")
df5 <- pivot_expr(blood[[5]],GSE ="GSE76705")

df <- rbind(df1,df2,df3,df4,df5)

ggplot(df,aes(x=GSE,y=expr, fill=Disease)) + geom_boxplot() + theme_classic()

dfvar <- df %>% 
  group_by(genes,Disease) %>%
  summarise(var=var(expr), mean=mean(expr)) 

ggplot(dfvar,aes(x=Disease,y=var, fill = Disease)) + 
  geom_boxplot() + 
  theme_classic()

b <- ggplot(dfvar,aes(x=Disease,y=var, fill = Disease)) + 
  geom_boxplot() + 
  theme_classic()

s <- ggplot(dfvar,aes(x=mean,y=var, group = Disease ,color = Disease, alpha= 0.7)) + 
  geom_point() + 
  geom_smooth(size = .5, color = "black",linetype="dashed", aes (fill = Disease)) +
  theme_classic()

ggarrange(b, s)




# dfmean <- df %>% 
#   group_by(genes, GSE,Disease) %>%
#   summarise(mean=mean(expr), sd=sd(expr)) %>%
#   group_by(genes, GSE) %>%
#   summarise(var = var(mean))
# 
# ggplot(dfmean,aes(x=GSE,y=log(var))) + geom_boxplot() + theme_classic()

