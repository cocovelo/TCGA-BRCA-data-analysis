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
`ge_data <- assay(data)`

This gene expression data matrix was then fed into the edgeR pipeline for RNA-sequencing
analysis. To do this, the matrix was converted to an edgeR object with `my_dgelist <- DGEList(ge_data)`.
Then `norm_dgelist <- normLibSizes(my_dgelist)` was carried out which estimates the
normalisation needed across the samples. Next, I used `filterByExpr` to identify genes
which were not expressed across the majority of samples `filtered_ge <- filterByExpr(norm_dgelist)`.
Using `table(filtered_ge)` it was possible to see that 42,329 genes were filtered at this
stage, leaving 18,331 for analysis. 

To remove the genes recommended for exclusion due to low expression:
`retained_ge <- norm_dgelist[filtered_ge, , keep.lib.size = FALSE]`
`nrow(retained_ge)` confirmed that the retained_ge data matrix has 18,331 rows
representing the 18,331 genes with expression data.

To carry out the next steps in the exploratory analysis of the data, I used the
DESeq2 package. First, I saved the unstranded count data from my data set:
`raw_counts <- assays(data)$unstranded`
Then I converted the data object from a SummarizedExperiment object to a DESeqDataSet
object:
`dds <- DESeqDataSetFromMatrix(countData = raw_counts,`
`colData = colData(data),`
`design = ~1)`
Next I performed a variance stabilising transformation to the data which is a normalisation
process:
`vst_data <- vst(dds, blind = TRUE)`
Before plotting the PCA plot: `plotPCA(vst_data, intgroup = "sample_type")`.
This plotted the samples on PC1 and PC2 and coloured the data points by sample type
i.e. whether the sample was metastatic, primary tumour or normal tissue.

![PCA Plot of TCGA-BRCA RNA-Seq Data](figures/pca_plot_tcga_brca.png)

## Clustering of samples




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

