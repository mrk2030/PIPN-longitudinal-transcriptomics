library(GEOquery)
getwd()
getGEOSuppFiles(
  "GSE185084",
  makeDirectory = FALSE,
  baseDir = "data_raw"
)
list.files("data_raw")