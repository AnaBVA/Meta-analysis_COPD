library(tidyverse)

tissue <- read_csv("Tissue/output_data/2020-10-07_Step3_Summary.csv")[-10,] #rm Argentina
blood <- read_csv("Blood/output_data/Blood_2020-11-25_Step3_Summary.csv")

tdf <- tissue %>% 
  select(GSE,CONTROL,COPD) %>%
  pivot_longer(!GSE, names_to = "disease", values_to = "count")

ggplot(tdf, aes(x=disease,y = count, fill = disease)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ GSE, ncol=4) + 
  ylab("Number of samples") +
  xlab("") +
  #scale_fill_manual(values=c("#66CC99","#9999CC")) +
  theme_bw()


 
bdf <- blood %>% 
  select(GSE,CONTROL,COPD) %>%
  pivot_longer(!GSE, names_to = "disease", values_to = "count")

ggplot(bdf, aes(x=disease,y = count, fill = disease)) +
  geom_bar(stat = "identity") +
  facet_wrap( ~ GSE, ncol=2) + 
  ylab("Number of samples") +
  xlab("") +
  #scale_fill_manual(values=c("#66CC99","#9999CC")) +
  theme_bw()


t <- tdf %>% 
  group_by(disease) %>%
  mutate(sum = sum(count),.keep = "none") %>%
  add_column(sample = "Lung tissue")  %>%
  unique()

b <- bdf %>% 
  group_by(disease) %>%
  mutate(sum = sum(count),.keep = "none") %>%
  add_column(sample = "Blood")  %>%
  unique()

df <- rbind(t,b)

ggplot(df, aes(x = sample, y = sum, fill = disease)) +
  geom_bar(position="dodge",stat = "identity")  +
  xlab("") +
  ylab("Total number of samples") +
  theme_bw(20)
