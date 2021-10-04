### Blood vs Lung tissue

library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)

#############
## meta analysis resulrs
blood_meta <- read_csv(here::here("Blood/output_data/Blood_2021-03-01_Step6_meta-analysis.csv"))
lung_meta <- read_csv(here::here("Tissue/output_data/2021-03-01_Step6_meta-analysis.csv"))

# filter genes with a qvalue
qval_cutoff <- 0.05
blood <- filter(blood_meta,qval.random <= qval_cutoff & I2 < 0.40 & num_exp > 3) %>% rename(genes = X1)
lung <- filter(lung_meta,qval.random <= qval_cutoff & I2 < 0.40 & num_exp >= 9) %>% rename(genes = X1)

#############
# Scatter plot

lfc <- left_join(lung_meta, blood_meta, by = "X1") %>% 
  #filter(I2.x < 0.4 & I2.y < 0.4) %>% 
  #filter(qval.random.x <= qval_cutoff & qval.random.y <= qval_cutoff) %>% 
  select(X1,TE.random.x, TE.random.y)

pdf(here::here("Lung-Blood/fig/scatter_lungvsblood.pdf"))
ggplot(lfc, aes(x = TE.random.x, y = TE.random.y)) +
  geom_point() +
  xlab("Lung tissue")+
  ylab("Blood") + 
  theme_classic(16)
dev.off()


#############

# data for the upset plot
us <- full_join(blood[,1:2], lung[,1:2], by= "genes") %>%
  mutate_if(is.numeric, ~1 * (. != "NA" )) %>% # change values to 1 if non NA
  replace(is.na(.), 0)

colnames(us) <- c("genes", "Blood", "Lung")

# upset plot
m1 = make_comb_mat(us)
ht <- draw(UpSet(m1,
                 pt_size = unit(3, "mm"), 
                 lwd = 3
))

od = column_order(ht)
cs = comb_size(m1)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(1, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 12))
})

#############
# feature selection for the model results
blood.rf_rfe <- readRDS("Blood/output_data/rf_rfe_701g.RDS")
blood.rf_rfe <- as_tibble(blood.rf_rfe$optVariables)
blood.rf_rfe$Blood <- 1

lung.rf_rfe <- readRDS("Tissue/output_data/rf_rfe_61g.RDS")
lung.rf_rfe <- as_tibble(lung.rf_rfe$optVariables)
lung.rf_rfe$Lung <- 1

# data for the upset plot
us <- full_join(blood.rf_rfe, lung.rf_rfe, by= "value") %>%
  mutate_if(is.numeric, ~1 * (. != "NA" )) %>% # change values to 1 if non NA
  replace(is.na(.), 0)

colnames(us) <- c("genes", "Blood", "Lung")

# upset plot
m1 = make_comb_mat(us)
ht <- draw(UpSet(m1,
                 pt_size = unit(3, "mm"), 
                 lwd = 3
))

od = column_order(ht)
cs = comb_size(m1)
decorate_annotation("intersection_size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(1, "pt"), 
            default.units = "native", just = "bottom", gp = gpar(fontsize = 12))
})

filter(us, Blood == Lung)

###########


## full tables with logFC
fblood <- read_csv(here::here("Blood/output_data/Blood_2020-11-28_Step5_Full_Tables.csv"))
flung <- read_csv(here::here("Tissue/output_data/2020-10-08_Step5_Full_Tables.csv"))

# all genes
bt <- filter(blood, blood$genes %in% lung$genes) %>% 
  select(genes)

# UP in blood and lung
up_b <- blood %>% 
  filter(TE.random > 0)

up_l <- lung %>% 
  filter(TE.random > 0)

up_bt <- filter(up_b, up_b$genes %in% up_l$genes) %>%
  select(genes)

# Immune response... oir is 25 shared genes in IR hallmark pathway
#bt <- bt[bt$genes %in% names(oir),]

# genes shared in blood and tissue
bt <- inner_join(blood, lung, by ="genes", suffix = c(".blood", ".lung"))
bt <- bt  %>%
  filter(tau2.blood != "NA" & tau2.lung != "NA") %>%
  rowwise() %>% mutate(TE.median = median(c(TE.random.blood,TE.random.lung))) %>%
  arrange(TE.median)

#Top genes
# bt <- rbind(head(bt,15), tail(bt,15))

#write.csv(bt, str_c("output_data/2021-04-10_shared_genes_q",qval_cutoff,".csv"), quote = F)

# filter genes in full table
fblood <- filter(fblood, fblood$GENE.SYMBOL %in% bt$genes) 
flung <- filter(flung, flung$GENE.SYMBOL %in% bt$genes)

## heatmap
# join tables
topgenes <- inner_join(fblood,flung, by= "GENE.SYMBOL")
topgenes <- topgenes %>%
  select("GENE.SYMBOL",str_subset(colnames(topgenes),"logFC")) %>%
  column_to_rownames("GENE.SYMBOL") 

topgenes <- topgenes[which(rowSums(is.na(topgenes)) < 9),]

topgenes <- scale(topgenes)

color <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(15))
breaks <- c(seq(min(topgenes,na.rm = T),-0.1,length.out = 7),0,seq(0.1,max(topgenes,na.rm = T),length.out = 7))

type <-  data.frame(type = c(rep("Blood",5),rep("Lung",11)))
row.names(type) <- c(str_subset(colnames(fblood),"logFC"),str_subset(colnames(flung),"logFC"))

my_colour = list(
  type = c(Blood = "#c45a75", Lung = "#5abdc4")
)

title <- str_c("Genes with q-value less than ",qval_cutoff, " (ngenes = ",dim(topgenes)[1], ")")

pheatmap(topgenes,
         main = title,
         color = color,
         breaks = breaks,
         #fontsize_row = 10,
         fontsize_col = 10,
         cluster_rows = T,
         clustering_distance_rows = "canberra",
         clustering_method = "ward.D2",
         cluster_cols = F,
         gaps_col = 5,
         cutree_rows =4,
         annotation_colors = my_colour,
         annotation_col = type,
         show_rownames = T,
         show_colnames = T,
         border_color = NA)


pheatmap(topgenes,
         main = "Top genes in blood and lung tissue",
         color = color,
         breaks = breaks,
         fontsize_row = 10,
         fontsize_col = 10,
         cluster_rows = T,
         clustering_distance_rows = "canberra",
         clustering_method = "ward.D",
         cluster_cols = F,
         gaps_col = 5,
         cutree_rows = 2,
         annotation_colors = my_colour,
         annotation_col = type,
         #show_rownames = F,
         show_colnames = T,
         border_color = NA)


# pheatmap(topgenes,
#          main = "Inflammatory Response genes in blood and lung tissue",
#          color = color,
#          breaks = breaks,
#          fontsize_row = 10,
#          fontsize_col = 10,
#          cluster_rows = T,
#          clustering_distance_rows = "manhattan",
#          clustering_method = "ward.D",
#          cluster_cols = F,
#          gaps_col = 5,
#          cutree_rows = 4,
#          annotation_colors = my_colour,
#          annotation_col = type,
#          #show_rownames = F,
#          show_colnames = T,
#          border_color = NA)
