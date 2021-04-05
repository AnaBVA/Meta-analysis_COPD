#########################################################################
## This is script is to find overlap genes in blood and lung using GTEx data v8
#########################################################################
## Download data
####################

## Annotation:
# wget https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt

## Transcriptomic data (Gene read counts)
# https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
#########################################################################
library(tidyr)
#########################################################################
# GTEx Data
gtex <-  read_tsv("Lung-Blood/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct",skip = 2)
a_gtex <- read_tsv("Lung-Blood/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")

###############################################
# Select Blood and Lung samples
l_sam <- a_gtex %>% 
  filter(a_gtex$SMTSD == "Lung" ) %>% 
  select(SAMPID) %>% 
  filter(SAMPID %in% colnames(gtex))

b_sam <- a_gtex %>% 
  filter(a_gtex$SMTSD == "Whole Blood" ) %>% 
  select(SAMPID) %>% 
  filter(SAMPID %in% colnames(gtex))

###############################################
# Select Blood and Lung gene expression data
l_gtex <- gtex %>% 
  select(c(Description,l_sam$SAMPID))

b_gtex <- gtex %>% 
  select(c(Description,b_sam$SAMPID))

# Remove complete gtex data frame
rm(gtex)

###############################################
# Summarizing values from same tissue samples per gene

l_median <- l_gtex %>% 
  rowwise(Description) %>%  # filtering gene names 
  mutate(m = median(c_across(where(is.numeric)))) %>%  # calculating the median per gene (row)
  select(Description,m)

b_median <- b_gtex %>% 
  rowwise(Description) %>%  # filtering gene names 
  mutate(m = median(c_across(where(is.numeric)))) %>% 
  select(Description,m)

###############################################
# remove genes without expression (!counts == 0)

l <- l_median %>% filter(m != 0) 
b <- b_median %>% filter(m != 0) 

###############################################
# Shared genes percentage between blood and lung

sum(l$Description %in% b$Description) / dim(l_gtex)[1]
# l %>% filter(Description %in% b$Description)


###############################################
## UpSet plot 

# data for the upset plot
s <- full_join(l,b, by= "Description",  suffix = c(".lung", ".blood")) %>%
  mutate_if(is.numeric, ~1 * (. != "NA" )) %>% # change values to 1 if non NA
  replace(is.na(.), 0)

colnames(s) <- c("genes", "Lung", "Blood")

# plot
m1 = make_comb_mat(s)
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



