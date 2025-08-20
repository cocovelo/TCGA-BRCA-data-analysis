## Exploratory data analysis
library(PCAtools)
library(edgeR)
library(DESeq2)
library(PCAtools)

saved_file_path <- "C:/Users/colin/Documents/R projects/RNA-seq-clustering/GDCdata/TCGA-BRCA-RNASeq-SummarizedExperiment.RData"
load(saved_file_path)
ge_data <- assay(data)

class(data) # should be class "SummarizedExperiment"
head(data)
dim(data)
raw_rownames <- rownames(se_object)
raw_colnames <- colnames(se_object)
grepl("ENSG", raw_rownames) %>% table()
grepl("TCGA", raw_colnames) %>% table()

saved_file_path <- "C:/Users/colin/Documents/R projects/RNA-seq-clustering/GDCdata/TCGA-BRCA-RNASeq-SummarizedExperiment.RData"
load(saved_file_path)
# This saves the se_object as the variable "data"

