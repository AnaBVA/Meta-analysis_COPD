library(tidyverse)
library(igraph)
library(RColorBrewer)
library(colourvalues)

# Data
install.packages("colourvalues")
df_gene <- read_delim(here::here("Lung-Blood/data/variant_gene_tf.csv"))

# Filter data 
df_fil <- df_gene %>% 
  filter(!(gene == "TP53" & abs(TF_TE.random) < 0.05)) %>% 
  filter(!(gene == "DSP" & abs(TF_TE.random) < 0.05))
  #filter(abs(TF_TE.random) > quantile(abs(df_gene$TF_TE.random))[2])
#df_fil <- df_gene

# Set graph vertices
df <- data.frame(from = df_fil$TF, to = df_fil$gene)

# Create relations
net <- graph_from_data_frame(d=df,  directed=T)
net <- simplify(net, remove.multiple = F, remove.loops = T)

# Vertice data
ver <- as_data_frame(net,what = "vertices")

# Is it a TF? 
V(net)$type <-  ver$name %in% df_gene$TF

# Shape
V(net)$shape <- c("circle","vrectangle")[V(net)$type+1]

#Color
#V(net)$color <- c("steel blue", "orange")[V(net)$type+1]
logfc <- data.frame(gene = c(df_fil$TF,df_fil$gene),
                    logfc = c(df_fil$TF_TE.random,df_fil$gene_TE.random))
logfc <- unique(logfc)
rownames(logfc) <- logfc$gene
V(net)$color <- colour_values(-logfc[ver$name,"logfc"], palette = "RdBu")

# Size
# deg <- degree(net, mode="all")
# V(net)$size <- log10(deg+1)
V(net)$label.cex = (7-V(net)$type) *.06

# Edge size
E(net)$arrow.size <- .15
#E(net)$edge.color <- "gray80"

# Topology
#layout_nicely
#layout_with_dh
l <- layout_nicely(net)
l <- norm_coords(l, ymin=-1, ymax=1, xmin=-1, xmax=1)

pdf(here::here("Lung-Blood/fig/var_graph.pdf"))
plot(net,
     vertex.label.color="black",
     vertex.size=(3-V(net)$type) *3.8, # Node's size 
     vertex.size2=(2-V(net)$type) *6,
     layout=l * 1,
     rescale = F,
     edge.arrow.mode=2,
     edge.color="gray40" #vertex.color="gray50"
     )
dev.off()
