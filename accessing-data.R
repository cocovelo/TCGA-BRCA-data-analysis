library(tidyverse)
library(GEOquery)
library(TCGAbiolinks)
library(SummarizedExperiment)

manifest <- read.delim("C:/Users/colin/Downloads/gdc_manifest.2025-08-18.140242.txt")

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = manifest$id
)

# Using the query object, then used the following in a terminal to download the files associated with gdc_manifest.txt:
# "C:/Users/colin/Downloads/gdc-client_2.3_Windows_x64-py3.8-windows-2019/gdc-client_2.3_Windows_x64/gdc-client.exe" download -m "C:/path/to/your/gdc_manifest.txt"
# This was necessary to download the data rather than using GDCdownload as there was a mismatch between
# the files included in my manifest and those expected by GDCdownload()

data_dir <- "C:/Users/colin/Documents/R projects/RNA-seq-clustering/GDCdata/"
output_filename <- "TCGA-BRCA-RNASeq-SummarizedExperiment.RData"

se_object <- GDCprepare(
  query = query,
  save = TRUE,
  save.filename = file.path(data_dir, output_filename),
  directory = data_dir
)

ge_data <- assay(se_object)
