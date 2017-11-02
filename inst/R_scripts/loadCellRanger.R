# Example script for Mus musculus
library(bcbioSingleCell)

# Cell Ranger reference datasets: https://goo.gl/3kqSu3
# Check the website for latest reference dataset version. These reference
# datasets are quite large and should be saved in a central location rather
# than per project.
dir.create("annotations", showWarnings = FALSE)
remoteFile <- file.path(
    "http://cf.10xgenomics.com",
    "supp",
    "cell-exp",
    "refdata-cellranger-mm10-1.2.0.tar.gz")
localFile <- file.path(
    "annotations",
    "refdata-cellranger-mm10-1.2.0.tar.gz")
download.file(remoteFile, localFile)
untar(localFile)
rm(remoteFile, localFile)

bcb <- loadCellRanger(
    uploadDir = file.path("data", "cellranger"),
    refDataDir = file.path("annotations", "refdata-cellranger-mm10-1.2.0"),
    sampleMetadataFile = file.path("meta", "sample_metadata.xlsx"),
    interestingGroups = c("genotype", "age")
)
# Back up all data inside bcbio object
flatFiles <- flatFiles(bcb)
saveData(bcb, flatFiles, dir = "data")
