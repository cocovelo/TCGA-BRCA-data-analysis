## Exploratory data analysis
library(PCAtools)
library(edgeR)
library(DESeq2)
library(PCAtools)
library(pheatmap)

saved_file_path <- "C:/Users/colin/Documents/R projects/RNA-seq-clustering/GDCdata/TCGA-BRCA-RNASeq-SummarizedExperiment.RData"
load(saved_file_path)
# This saves the se_object as the variable "data"
raw_counts_matrix <- assays(data)$unstranded
sample_data_df <- as.data.frame(colData(data))

class(data) # should be class "SummarizedExperiment"
head(data)
dim(data)

dds <- DESeqDataSetFromMatrix(countData = raw_counts_matrix,
                              colData = sample_data_df,
                              design = ~1)

keep <- filterByExpr(dds)
dds_filtered <- dds[keep,]
print(paste("Original number of genes:", nrow(dds)))
print(paste("Number of genes after filtering:", nrow(dds_filtered)))
vst_data_filtered <- vst(dds_filtered, blind = TRUE)
plotPCA(vst_data, intgroup = "sample_type")

sampleDists <- dist(t(assay(vst_data_filtered)))
sampleDistMatrix <- as.matrix(sampleDists)
annotation_col <- as.data.frame(colData(vst_data_filtered)[, c("sample_type", "gender")])
pheatmap(sampleDistMatrix,
  clustering_distance_rows=sampleDists,
  clustering_distance_cols=sampleDists,
  annotation_col=annotation_col,
  main="Sample-to-Sample Distance Heatmap",
  show_rownames = FALSE,
  show_colnames = FALSE
)
