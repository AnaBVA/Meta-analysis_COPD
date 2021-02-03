### Blood vs Lung tissue

library(tidyverse)
library(ComplexHeatmap)
library(pheatmap)
library(RColorBrewer)

## meta analysis resulrs
blood <- read_csv(here::here("Blood/output_data/Blood_2021-02-01_Step6_meta-analysis.csv"))
lung <- read_csv(here::here("Tissue/output_data/2021-02-02_Step6_meta-analysis.csv"))

# filter genes with a qvalue
qval_cutoff <- 0.05
blood <- filter(blood,qval.fixed <= qval_cutoff) %>% rename(genes = X1)
lung <- filter(lung,qval.fixed <= qval_cutoff) %>% rename(genes = X1)


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


## full tables with logFC
fblood <- read_csv(here::here("Blood/output_data/Blood_2020-11-28_Step5_Full_Tables.csv"))
flung <- read_csv(here::here("Tissue/output_data/2020-10-08_Step5_Full_Tables.csv"))

# genes shared in blood and tissue
bt <- inner_join(blood, lung, by ="genes", suffix = c(".blood", ".lung"))
bt <- bt  %>%
  filter(tau2.blood != "NA" & tau2.lung != "NA") %>%
  rowwise() %>% mutate(TE.median = median(c(TE.fixed.blood,TE.fixed.lung))) %>%
  arrange(TE.median)

 bt <- rbind(head(bt,15), tail(bt,15))

#bt <- filter(blood, blood$genes %in% lung$genes) %>% 
#  select(genes)

# filter genes in full table
fblood <- filter(fblood, fblood$GENE.SYMBOL %in% bt$genes) 
flung <- filter(flung, flung$GENE.SYMBOL %in% bt$genes)

## heatmap
# join tables
topgenes <- inner_join(fblood,flung, by= "GENE.SYMBOL")
topgenes <- topgenes %>%
  select("GENE.SYMBOL",str_subset(colnames(topgenes),"logFC")) %>%
  column_to_rownames("GENE.SYMBOL") 

topgenes <- topgenes[which(rowSums(is.na(topgenes)) < 11),]

topgenes <- scale(topgenes)

color <- rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(15))
breaks <- c(seq(min(topgenes,na.rm = T),-0.1,length.out = 7),0,seq(0.1,max(topgenes,na.rm = T),length.out = 7))

type <-  data.frame(type = c(rep("Blood",5),rep("Lung",11)))
row.names(type) <- c(str_subset(colnames(fblood),"logFC"),str_subset(colnames(flung),"logFC"))

my_colour = list(
  type = c(Blood = "#851d1d", Lung = "#1d852e")
)

pheatmap(topgenes,
         main = "Genes in blood and lung tissue",
         color = color,
         breaks = breaks,
         #fontsize_row = 10,
         fontsize_col = 10,
         cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D",
         cluster_cols = F,
         gaps_col = 5,
         cutree_rows = 4,
         annotation_colors = my_colour,
         annotation_col = type,
         show_rownames = F,
         show_colnames = T,
         border_color = NA)


pheatmap(topgenes,
         main = "Top genes in blood and lung tissue",
         color = color,
         breaks = breaks,
         fontsize_row = 10,
         fontsize_col = 10,
         cluster_rows = T,
         clustering_distance_rows = "euclidean",
         clustering_method = "ward.D",
         cluster_cols = F,
         gaps_col = 5,
         cutree_rows = 2,
         annotation_colors = my_colour,
         annotation_col = type,
         #show_rownames = F,
         show_colnames = T,
         border_color = NA)
