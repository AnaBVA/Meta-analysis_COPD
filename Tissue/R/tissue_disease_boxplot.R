library(tidyverse)

# Data
tissue <- readRDS(here::here("Tissue/output_data/2020-10-07_Step3_LungTissue-CURATED.RDS"))
#fc_tissue <- read_csv("Tissue/output_data/2020-10-07_Step5_Full_Tables.csv")

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
  df$expr <- as.numeric(df$expr)
  
  return(df)
}

df1 <- lapply(c(1:length(tissue))[-c(4,10,12)], # 4 and 12 are RNAseq, 10 is low quality
             function(i){pivot_expr(tissue[[i]][[1]],GSE = names(tissue[i]))}
             )

df <- bind_rows(df1)

ggplot(df,aes(x=GSE,y=expr, fill=Disease)) + 
  geom_boxplot() + 
  #coord_cartesian(ylim = c(0,30)) +
  theme_classic() 
