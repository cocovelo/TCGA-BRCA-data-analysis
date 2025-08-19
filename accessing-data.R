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
# the files included in my manifest and those expected b
se_object <- GDCprepare(
  query = query,
  save = TRUE,
  )

