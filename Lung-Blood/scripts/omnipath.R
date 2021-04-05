## This script is to use omnipathR for createing a protein-protein interaction

library(OmnipathR)
library(tidyverse)
library(dnet)
library(gprofiler2)

# meta-analysis results from blood
bmeta <- read_csv(here::here("Blood/output_data/Blood_2021-03-01_Step6_meta-analysis.csv"))

# meta-analysis results from lung
lmeta <- read_csv(here::here("Tissue/output_data/2021-03-01_Step6_meta-analysis.csv"))


# select genes with a qval less than 0.05
bmeta <- bmeta %>% 
  filter(qval.random <= 0.05 & num_exp > 3 & I2 < .40) %>% 
  arrange(desc(TE.random))

lmeta <- lmeta %>% 
  filter(qval.random <= 0.05 & num_exp >= 9 & I2 < .40) %>% 
  arrange(desc(TE.random))


# all genes
bt <- filter(bmeta, bmeta$X1 %in% lmeta$X1) %>% 
  select(X1) 


## The interactions are stored into a data frame.
## The interactions are stored into a data frame.
interactions <-
  import_omnipath_interactions(resources=c("connectomeDB2020",
                                           "STRING_ICELLNET", "STRING_talklr", 
                                           "Reactome_ICELLNET", "Reactome_ProtMapper", "Reactome_SignaLink3",
                                           "DoRothEA", "RegNetwork_DoRothEA"))


## We select all the interactions in which Amfr gene is involved
interactions_Amfr <- filter(interactions, source_genesymbol %in% bt$X1|
                                     target_genesymbol %in% bt$X1)

## We print these interactions:
print_interactions(interactions_Amfr)

## We transform the interactions data frame into a graph
OPI_g <- interaction_graph(interactions = interactions_Amfr)

OPI_g_undirected <- as.undirected(OPI_g, mode=c("mutual"))
OPI_g_undirected <- simplify(OPI_g_undirected)
cl_results <- cluster_fast_greedy(OPI_g_undirected)
## We extract the cluster where a protein of interest is contained
cluster_id <- cl_results$membership
module_graph <- induced_subgraph(OPI_g,
                                 V(OPI_g)$name[which(cl_results$membership == cluster_id)])


## We print that cluster with its interactions.
# par(mar=c(0.1,0.1,0.1,0.1))
plot(module_graph, vertex.label.color="black",vertex.frame.color="#ffffff",
     vertex.size= 15, edge.curved=.2,
     vertex.color = ifelse(igraph::V(module_graph)$name == "ERBB2","yellow",
                           "#00CCFF"), edge.color="blue",edge.width=0.8)


