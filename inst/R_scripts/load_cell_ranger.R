# nolint start
#
# Load 10X Genomics Cell Ranger Data
#
# Michael Steinbaugh
# 2018-02-26
#
# Latest version of this script is available here:
# script <- system.file(
#     "R_scripts/load_cell_ranger.R",
#     package = "bcbioSingleCell")
# file.edit(script)
#
# nolint end

library(bcbioSingleCell)

# Cell Ranger reference datasets: https://goo.gl/3kqSu3
# Check the website for latest reference dataset version. These reference
# datasets are quite large and should be saved in a central location rather
# than per project.
dir_create("annotations")
remote <- path(
    "http://cf.10xgenomics.com",
    "supp",
    "cell-exp",
    "refdata-cellranger-mm10-1.2.0.tar.gz")
local <- path(
    "annotations",
    "refdata-cellranger-mm10-1.2.0.tar.gz")
download.file(remote, local)
untar(local)
rm(remote, local)

bcb <- loadCellRanger(
    uploadDir = path("data", "cellranger"),
    refDataDir = path("annotations", "refdata-cellranger-mm10-1.2.0"),
    sampleMetadataFile = path("meta", "sample_metadata.xlsx"),
    interestingGroups = c("genotype", "age")
)
saveData(bcb, dir = path("data", Sys.Date()))
