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
the files included in my manifest and those expected by `GDCdownload()`. I then prepared the data in R using:

`se_object <- GDCprepare(
        query = query,
        save = TRUE,
        )`

