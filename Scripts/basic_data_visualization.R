# library load


library(tidyverse)
library(ggplot2)

# load data

data_long <- read.csv("data/GSE183947_long.csv")


# bar plot

data_long |> 
  filter(gene=="BRCA1") |> 
  ggplot(aes(x= sample, y=FPKM, fill = Tissue)) +
  geom_col()
  
# density plot

data_long |> 
  filter(gene=="BRCA1") |> 
  ggplot(aes(x=FPKM, fill = Tissue)) +
  geom_density(alpha = 0.3)
 

# box plot

data_long |> 
  filter(gene=="BRCA1") |> 
  ggplot(aes(x=Metastasis,y=FPKM,fill = Tissue )) +
  geom_boxplot()

# violin plot

data_long |> 
  filter(gene=="BRCA1") |> 
  ggplot(aes(x=Metastasis,y=FPKM )) +
  geom_violin()


# scatter plot

data_long |> 
  filter(gene == "BRCA1"| gene=="BRCA2") |> 
  spread(key=gene, value = FPKM) |> 
  ggplot(aes(x=BRCA1,y=BRCA2, colour = Tissue))+
  geom_point()+
  geom_smooth(method = "lm", se = F)
 

# heatmap


gene_of_interest <- c("BRCA1", "BRCA2", "TP53", "ALK", "MYCN")

data_long |> 
  filter(gene %in% gene_of_interest) |> 
  ggplot(aes(x=sample, y= gene, fill = FPKM))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")


# the darker the color the lower the expression
# the lower the color the greater the expression


# saving data with command line

gene_of_interest <- c("BRCA1", "BRCA2", "TP53", "ALK", "MYCN")

p <- data_long |> 
  filter(gene %in% gene_of_interest) |> 
  ggplot(aes(x=sample, y= gene, fill = FPKM))+
  geom_tile()+
  scale_fill_gradient(low = "white", high = "red")

ggsave(p, filename = "figure/prr.pdf", width = 10, height = 9)























