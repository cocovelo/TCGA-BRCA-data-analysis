# Accessing and downloading RNA-seq data from GEO181466 and GEO135298
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181466
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135298

library(tidyverse)
library(GEOquery)
library(TCGAbiolinks)
library(SummarizedExperiment)

file_path_gse181466 <- "C:/Users/colin/Downloads/GSE181466_series_matrix.txt/GSE181466_series_matrix.txt"
file_path_gse135298 <- "C:/Users/colin/Downloads/GSE135298_series_matrix.txt/GSE135298_series_matrix.txt"

all_lines_gse181466 <- read_lines(file_path_gse181466)
all_lines_gse135298 <- read_lines(file_path_gse135298)
metadata_gse181466 <- all_lines_gse181466[grep("^!", all_lines_gse181466)]
metadata_gse135298 <- all_lines_gse135298[grep("^!", all_lines_gse135298)]
head(metadata_gse181466)
head(metadata_gse135298)

## reading in gse181466 PAM50 subtype information

characteristics_line <- all_lines_gse181466[grep("!Sample_characteristics_ch1", all_lines_gse181466)]
subtype_line <- characteristics_line[grep("subtype:", characteristics_line)]
subtype_data <- strsplit(subtype_line, "\t")[[1]]

clean_subtypes <- gsub("subtype: ", "", subtype_data)
clean_subtypes <- gsub("\"", "", clean_subtypes)

data_header <- all_lines_gse181466[grep("!Sample_geo_accession", all_lines_gse181466)]
sample_ids_181466 <- strsplit(data_header, "\t")[[1]]
sample_ids_181466 <- sample_ids_181466[-1]
sample_ids_181466 <- gsub("\"", "", sample_ids_181466)

pam50_metadata_181466 <- data.frame(
SampleID = sample_ids_181466,
PAM50_Subtype = clean_subtypes[-1]
)

## reading in gse135298 PAM50 subtype information

characteristics_line_135298 <- all_lines_gse135298[grep("!Sample_characteristics_ch1", all_lines_gse135298)]
pam50_line <- characteristics_line_135298[grep("pam50:", characteristics_line_135298)]
pam50_data <- strsplit(pam50_line, "\t")[[1]]
clean_pam50 <- gsub("pam50: ", "", pam50_data)
clean_pam50 <- gsub("\"", "", clean_pam50)
geo_accession_line <- all_lines_gse135298[grep("!Sample_geo_accession", all_lines_gse135298)]
sample_ids_135298 <- strsplit(geo_accession_line, "\t")[[1]]
sample_ids_135298 <- sample_ids_135298[-1]
sample_ids_135298 <- gsub("\"", "", sample_ids_135298)
pam50_metadata_135298 <- data.frame(
  SampleID = sample_ids_135298,
  PAM50_Subtype = clean_pam50[-1]
)

## reading in gse181466 gene expression data

gse181466_expression <- read_delim(
  "C:/Users/colin/Downloads/GSE181466_rsem_genes_norm_matrix-97.txt/GSE181466_rsem_genes_norm_matrix-97.txt", 
  delim = "\t", 
  skip = 89, 
  col_names = FALSE
)
gse181466_expression <- as.data.frame(gse181466_expression)
rownames(gse181466_expression) <- gse181466_expression[, 1]
gse181466_expression <- gse181466_expression[, -1]
colnames(gse181466_expression) <- sample_ids_181466
head(gse181466_expression)

## reading in gse135298 gene expression data

data_start_line <- grep("^\"ID_REF\"", all_lines_gse135298)
gse135298_expression <- read_delim(
  "C:/Users/colin/Downloads/GSE135298_OSLO2_EMIT0_RNA-seq.txt/GSE135298_OSLO2_EMIT0_RNA-seq.txt",
  delim = "\t",
  skip = data_start_line,
  col_names = FALSE
)
gse135298_expression <- as.data.frame(gse135298_expression)
rownames(gse135298_expression) <- gse135298_expression$ID_REF
gse135298_expression$ID_REF <- NULL
head(gse135298_expression)
rownames(gse135298_expression) <- gse135298_expression$X1
gse135298_expression$X1 <- NULL
colnames(gse135298_expression) <- sample_ids_135298
head(gse135298_expression)


## next steps are to
## perform vst-normalisation on the data sets
## combine the metadata and normalised gene expression tables for each gse study
## then perform the machine learning steps

table(pam50_metadata_135298$SampleID == colnames(gse135298_expression))
table(pam50_metadata_181466$SampleID == colnames(gse181466_expression))

