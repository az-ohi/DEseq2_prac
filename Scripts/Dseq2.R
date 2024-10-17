library(DESeq2)
library(tidyverse)
library(airway)

# read count data

counts_data <- read.csv("data/counts_data.csv")
head(counts_data)

# read in sample info

codData <- read.csv("data/sample_info.csv")

# checking

all(colnames(counts_data) %in% rownames(codData))

# same order?

all(colnames(counts_data) == rownames(codData))


# construct deseqDataSet object

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = codData,
                       design = ~ dexamethasone
                       )
dds

# filtering low count gene

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

dds









