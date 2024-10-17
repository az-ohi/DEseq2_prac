library(DESeq2)
library(tidyverse)
library(airway)

# read count data

counts_data <- read.csv("data/counts_data.csv")
head(counts_data)

# read in sample info

colData <- read.csv("data/sample_info.csv")

# checking

all(colnames(counts_data) %in% rownames(colData))

# same order?

all(colnames(counts_data) == rownames(colData))


# construct deseqDataSet object

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~ dexamethasone
                       )
dds

# filtering low count gene

keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

dds


# set the factor level

dds$dexamethasone <- relevel(dds$dexamethasone, ref = "untreated" )

dds$dexamethasone


# NOTE!!! need to collapse technical replicates before run DEseq2
# run DEseq

dds <- DESeq(dds)





































