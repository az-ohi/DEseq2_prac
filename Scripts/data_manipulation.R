# library load

library(dplyr)
library(tidyverse)
library(GEOquery)

# read data

count_data <- read.csv("data/GSE183947_fpkm.csv")

count_data |> 
  head()

# get metadata

gse <- getGEO(GEO = "GSE183947", GSEMatrix = T)
metadata <- pData(phenoData(gse[[1]]))

# subsetting metadata

metadata_subset <- metadata |> 
  select(c(1,10,11,17))

# processing data

metadata_modified <- metadata_subset |> 
  rename(Tissue = characteristics_ch1, Metastasis = characteristics_ch1.1) |> 
  mutate(Tissue = gsub("tissue: ", "", Tissue)) |> 
  mutate(Metastasis = gsub("metastasis: ", "", Metastasis))


head(metadata_modified)

# reshaping data

data_long <- count_data |> 
  rename(gene = X) |> 
  gather(-gene,
               key = "sample",
               value = "FPKM"
               )

# join dataframe data.long+ metadata_modified

data_long <- data_long |> 
  left_join(metadata_modified, by=c("sample"="description"))


# export data

write.csv(data_long,"data/GSE183947_long.csv", row.names = F)


# explore data

data_long |> 
  filter(gene=="BRCA1"| gene=="BRCA2") |> 
  group_by(gene,Tissue) |>
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) |>
  arrange(desc(mean_FPKM))
  head()
    








































