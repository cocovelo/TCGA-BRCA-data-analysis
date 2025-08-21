# Identification, access, harvesting and analysis of a TCGA BRCA data set

## Identification of samples, downloading associated data and loading into R

Initially, I used the GDC data portal to filter samples and identify samples for use. Since I did this manually,
I downloaded the following from the GDC portal:

cohort.2025-08-18.tsv
gdc_sample_sheet.2025-08-18.tsv
files-table.2025-08-18.tsv
gdc-client_2.3_Windows_x64-py3.8-windows-2019
gdc_manifest.2025-08-18.140242
metadata.cart.2025-08-18.json
gdc_sample_sheet.2025-08-18(1).tsv
clinical.cart.2025-08-18.tar
biospecimen.cart.2025-08-18.tar

I next created a `manifest` R object with:
`manifest <- read.delim("C:/Users/colin/Downloads/gdc_manifest.2025-08-18.140242.txt")`

The next step was to create a query variable which describes the TCGA project:
`query <- GDCquery(
        project = "TCGA-BRCA",
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        barcode = manifest$id
        )`

Using the query object, I then used the following in a terminal to download the files associated with gdc_manifest.txt:
`"C:/Users/colin/Downloads/gdc-client_2.3_Windows_x64-py3.8-windows-2019/gdc-client_2.3_Windows_x64/gdc-client.exe" download -m "C:/path/to/your/gdc_manifest.txt"`

It was necessary to download the data in this way rather than using `GDCdownload()` as there was a mismatch between
the files included in my manifest and those expected by `GDCdownload()`. I then prepared the data in R using
GDCprepare, however in order to do this I had to first specify where the data were saved and what the
output file would be named:

`output_filename <- "TCGA-BRCA-RNASeq-SummarizedExperiment.RData"
data_dir <- "C:/Users/colin/Documents/R projects/RNA-seq-clustering/GDCdata/"`

I then created the `se_object` with:

`se_object <- GDCprepare(
    query = query,
    save = TRUE,
    save.filename = file.path(data_dir, output_filename),
    directory = data_dir
)`


## Exploratory data analysis

Now that the data have been successfully loaded into an R object, the next step is to
investigate the object to determine how many dimensions, what data are included etc.
To begin with, I used `dim(se_object)` which revealed 60,660 rows and 885 columns.
`str(se_object)` was impractical in this situation as there was too much information
to see in the console. Using `colnames(se_object)` and `rownames(se_object)` I was
able to discern that rows contained the gene expression data and other variables,
while the columns contained the samples. 

Sample names were, for example, in the format:
`colnames(se_object[1:10])`
  [1] "TCGA-GM-A2DL-01A-11R-A18M-07" "TCGA-AC-A2QI-01A-12R-A19W-07" "TCGA-EW-A1PD-01A-11R-A144-07"

While row names were, for example:
`rownames(se_object[1:10])`
 [1] "ENSG00000000003.15" "ENSG00000000005.6"  "ENSG00000000419.13" "ENSG00000000457.14" "ENSG00000000460.17"
 [6] "ENSG00000000938.13" "ENSG00000000971.16" "ENSG00000001036.14" "ENSG00000001084.13" "ENSG00000001167.14"


Row names and column names were saved:
`raw_rownames <- rownames(se_object)`
`raw_colnames <- colnames(se_object)`

I then attempted to determine how many of the rows contained gene expression data by
counting the number of times "ENSG" appeared in the name of a row:
`grepl("ENSG", raw_rownames) %>% table()`

This revealed that all 60,660 of the rows pertained to gene expression data, while
`grepl("TCGA", raw_colnames) %>% table()` revealed that all columns contained the
samples and no other variables.

`head(rowData(se_object))`, `head(colData(se_object))` and `class(se_object)` confirmed
that the se_object contained the additional metadata variables.

After closing the R session and re-opening a new session, I re-loaded the data into R
using

`saved_file_path <- "C:/Users/colin/Documents/R projects/RNA-seq-clustering/GDCdata/TCGA-BRCA-RNASeq-SummarizedExperiment.RData"`
`load(saved_file_path)`

Then confirmed the class was a "SummarizedExperiment" object using `class(data)`.
Then saved only the gene expression data as a new object:
`raw_counts_matrix <- assays(data)$unstranded`and the metadata/additional variables
`sample_data_df <- as.data.frame(colData(data))`.

The raw_counts_matrix was then converted to a DESeqDataSet object for processing of the
data:
`dds <- DESeqDataSetFromMatrix(countData = raw_counts_matrix,`
                              `colData = sample_data_df,`
                              `design = ~1)`

Next, I used `filterByExpr` to identify genes which were not expressed across the majority
of samples:
`keep <- filterByExpr(dds)` and `dds_filtered <- dds[keep,]`.
`print(paste("Original number of genes:", nrow(dds)))` and
`print(paste("Number of genes after filtering:", nrow(dds_filtered)))` confirmed that
there were 60,660 genes initially in the matrix but after filtering that number was
reduced to 18,303.

Next I performed a variance stabilising transformation to the data which is a normalisation
technique suitable for machine learning processes:
`vst_data_filtered <- vst(dds_filtered, blind = TRUE)`
Before plotting a PCA plot to visualise the samples and their relationships to one another:
`plotPCA(vst_data_filtered, intgroup = "sample_type")`.
This plotted the samples on PC1 and PC2 and coloured the data points by sample type
i.e. whether the sample was metastatic, primary tumour or normal tissue. Based on this PCA
plot there were no major outliers in the data set.

![PCA Plot of TCGA-BRCA RNA-Seq Data](figures/pca_plot_tcga_brca.png)

Next I investigated whether the type of cancer was associated with variation in the data
by plotting the PCA plot again, but this time grouped by "paper_BRCA_Subtype_PAM50":

`plotPCA(vst_data, intgroup = "paper_BRCA_Subtype_PAM50")`

![PCA plot grouped by PAM50](figures/pca-pam50.png)

This plot clearly showed that the subtype of cancer was a major source of variation in
the data. Another plot of PC1 score vs subtype confirmed this observation:

`ggplot(pcaData, aes(x = paper_BRCA_Subtype_PAM50, y = PC1, fill = paper_BRCA_Subtype_PAM50)) +`
`  geom_boxplot() +`
`  labs(title = "PC1 Scores by PAM50", x = "PAM50", y = "PC1 Score") +`
`  theme_minimal()`

![PC1 and PAM50](figures/pc1-pam50.png)

As an additional check, I wanted to ensure that the principal components were normally
distributed so that I could discern whether to use ANOVA or Kruskal-Wallis statistical
tests when investigating relationships between the PCs and other metadata variables.

`hist(pcaData$PC1)` and `hist(pcaData$PC2)` confirmed a normal distribution of both PC1 and PC2.

I then performed an ANOVA test between "sample-type" and PC1, and "sample-type" and PC2:
`aov_result_pc1 <- aov(PC1 ~ sample_type, data = pcaData)`
`summary(aov_result_pc1)`
`aov_result_pc2 <- aov(PC2 ~ sample_type, data = pcaData)`
`summary(aov_result_pc2)`
This confirmed that sample-type was statistically significantly associated with PC1 and PC2.

I performed the same tests for the "paper_BRCA_Subtype_PAM50" variable:
`aov_result_pc1 <- aov(PC1 ~ paper_BRCA_Subtype_PAM50, data = pcaData)`
`summary(aov_result_pc1)`
`aov_result_pc2 <- aov(PC2 ~ paper_BRCA_Subtype_PAM50, data = pcaData)`
`summary(aov_result_pc2)`
These tests confirmed that paper_BRCA_Subtype_PAM50  was also significantly associated with
PC1 and PC2. Testing of race confirmed that race was significantly associated with the PCs
which means that race should be considered in any future statistical modelling.

For continuous variables, I used linear regression to determine whether or not age at diagnosis
was significantly associated with PC1. Plotting the results of this suggested a positive
correlation between PC1 score and age at diagnosis, though the correlation was weak
(Pearson correlation = 0.177):

`ggplot(pcaData, aes(x = age_at_diagnosis, y = PC1)) +`
`    geom_point(alpha = 0.6) +`
`    geom_smooth(method = "lm", se = FALSE, color = "blue") + # Add a linear regression line`
`    labs(title = "PC1 Scores vs. Age at Diagnosis", x = "Age at Diagnosis", y = "PC1 Score") +`
`    theme_minimal()`

`if ("age_at_diagnosis" %in% colnames(pcaData)) {`
`  cor_pc1_age <- cor(pcaData$PC1, pcaData$age_at_diagnosis, use = "pairwise.complete.obs")`
`  print(paste("Pearson correlation (PC1 vs Age):", round(cor_pc1_age, 3)))`
`}`

![PC1 vs age at diagnosis scatter plot](figures/age-pc1-scatter-plot.png)

## Clustering of samples

Next I evaluated the relationships between the samples using a heatmap. To do this I used the
pheatmap package. I calculated the distances between the samples first using:
`sampleDists <- dist(t(assay(vst_data_filtered)))`
`sampleDistMatrix <- as.matrix(sampleDists)`
Then collected the metadata for annotation:
`annotation_col <- as.data.frame(colData(vst_data_filtered)[, c("sample_type", "gender")])`

I used the following to create the heatmap plot:
`pheatmap(sampleDistMatrix,`
`  clustering_distance_rows=sampleDists,`
`  clustering_distance_cols=sampleDists,`
`  annotation_col=annotation_col,`
`  main="Sample-to-Sample Distance Heatmap",`
`  show_rownames = FALSE,`
`  show_colnames = FALSE`
`)`

![Sample-to-sample distance heatmap](figures/sampletosampledistanceheatmap.png)

## Machine learning applications with pytorch and sklearn



## Deep learning



## Dependencies
library(tidyverse)
library(GEOquery)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(PCAtools)
library(edgeR)
library(DESeq2)
library(PCAtools)

