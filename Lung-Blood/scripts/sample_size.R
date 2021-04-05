library(tidyverse)
library(ggpubr)

tissue <- read_csv("Tissue/output_data/2020-10-07_Step3_Summary.csv")[-10,] #rm Argentina
blood <- read_csv("Blood/output_data/Blood_2020-11-25_Step3_Summary.csv")

tdf <- tissue %>% 
  select(GSE,CONTROL,COPD) %>%
  pivot_longer(!GSE, names_to = "disease", values_to = "count")

tp <- ggplot(tdf, aes(x=disease,y = count, fill = disease)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ GSE, ncol=4) + 
  ylab("Number of samples") +
  xlab("") +
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  geom_text(aes(label=count), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw(15)


 
bdf <- blood %>% 
  select(GSE,CONTROL,COPD) %>%
  pivot_longer(!GSE, names_to = "disease", values_to = "count")

bp <- ggplot(bdf, aes(x=disease,y = count, fill = disease)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ GSE, ncol=2) + 
  ylab("Number of samples") +
  xlab("") +
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  geom_text(aes(label=count), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw(15)


t <-   tdf %>% 
  group_by(disease) %>%
  mutate(sum = sum(count),.keep = "none",
         exp = length(unique(GSE))) %>%
  add_column(sample = "Lung tissue")  %>%
  unique()

b <- bdf %>% 
  group_by(disease) %>%
  mutate(sum = sum(count),.keep = "none",
         exp = length(unique(GSE))) %>%
  add_column(sample = "Blood")  %>%
  unique()

df <- rbind(t,b)

all <- ggplot(df, aes(x = sample, y = sum, fill = disease)) +
  geom_bar(position="dodge",stat = "identity")  +
  xlab("") +
  ylab("Total number of samples") +
  scale_fill_manual(values=c("#66CC99","#9999CC")) +
  geom_text(aes(label=sum), position=position_dodge(width=0.9), vjust=-0.25) +
  theme_bw(15)

d <- df %>% 
  group_by(sample) %>% 
  mutate(Samples = sum(sum),
         Experiments = as.double(exp)) %>% 
  select(sample,Experiments,Samples) %>% 
  unique() %>% 
  ungroup() %>% 
  add_row(sample = "Total",
          Experiments = sum(.$Experiments),
          Samples = sum(.$Samples)) %>% 
  rename(Sample_Type = sample)


dtabla <- ggtexttable(d) 

all + annotation_custom(ggplotGrob(dtabla),
                        xmin = 7)


ggpubr::ggarrange(tp,bp,all, dtabla,
                  labels = c("b", "c", "d"),
                  common.legend = TRUE,
                  heights = c(1.5,1),
                  widths = c(1.7,1),
                  nrow = 2, 
                  ncol = 2)

