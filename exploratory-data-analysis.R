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
plotPCA(vst_data_filtered, intgroup = "sample_type")
plotPCA(vst_data_filtered, intgroup = "paper_BRCA_Subtype_PAM50")
plotPCA(vst_data_filtered, intgroup = "ajcc_pathologic_stage")

pcaData <- plotPCA(vst_data_filtered, intgroup = "sample_type", returnData = TRUE)


ggplot(pcaData, aes(x = sample_type, y = PC1, fill = sample_type)) +
geom_boxplot() +
labs(title = "PC1 Scores by Sample Type", x = "Sample Type", y = "PC1 Score") +
theme_minimal()

ggplot(pcaData, aes(x = paper_BRCA_Subtype_PAM50, y = PC1, fill = paper_BRCA_Subtype_PAM50)) +
  geom_boxplot() +
  labs(title = "PC1 Scores by PAM50", x = "PAM50", y = "PC1 Score") +
  theme_minimal()

hist(pcaData$PC1)
hist(pcaData$PC2)
aov_result_pc1 <- aov(PC1 ~ sample_type, data = pcaData)
summary(aov_result_pc1)
aov_result_pc2 <- aov(PC2 ~ sample_type, data = pcaData)
summary(aov_result_pc2)

aov_result_pc1 <- aov(PC1 ~ paper_BRCA_Subtype_PAM50, data = pcaData)
summary(aov_result_pc1)
aov_result_pc2 <- aov(PC2 ~ paper_BRCA_Subtype_PAM50, data = pcaData)
summary(aov_result_pc2)

ggplot(pcaData, aes(x = age_at_diagnosis, y = PC1)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = "PC1 Scores vs. Age at Diagnosis", x = "Age at Diagnosis", y = "PC1 Score") +
    theme_minimal()

if ("age_at_diagnosis" %in% colnames(pcaData)) {
  cor_pc1_age <- cor(pcaData$PC1, pcaData$age_at_diagnosis, use = "pairwise.complete.obs")
  print(paste("Pearson correlation (PC1 vs Age):", round(cor_pc1_age, 3)))
}

## clustering
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

