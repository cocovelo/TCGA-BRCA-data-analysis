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
plotPCA(vst_data, intgroup = "paper_BRCA_Subtype_PAM50")

pcaData <- plotPCA(vst_data_filtered, intgroup = "sample_type", returnData = TRUE)


ggplot(pcaData, aes(x = sample_type, y = PC1, fill = sample_type)) +
geom_boxplot() +
labs(title = "PC1 Scores by Sample Type", x = "Sample Type", y = "PC1 Score") +
theme_minimal()

ggplot(pcaData, aes(x = paper_BRCA_Subtype_PAM50, y = PC1, fill = paper_BRCA_Subtype_PAM50)) +
  geom_boxplot() +
  labs(title = "PC1 Scores by PAM50", x = "PAM50", y = "PC1 Score") +
  theme_minimal()


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

