# nolint start
#
# Load 10X Genomics Cell Ranger Data
#
# Michael Steinbaugh
# 2018-03-12
#
# Latest version of this script is available here:
# script <- system.file(
#     "R_scripts/load_cell_ranger.R",
#     package = "bcbioSingleCell"
# )
# file.edit(script)
#
# nolint end

library(bcbioSingleCell)

# Cell Ranger reference datasets: https://goo.gl/3kqSu3
# Check the website for latest reference dataset version. These reference
# datasets are quite large and should be saved in a central location rather
# than per project.
dir.create("annotations", showWarnings = FALSE)
remote <- file.path(
    "http://cf.10xgenomics.com",
    "supp",
    "cell-exp",
    "refdata-cellranger-mm10-1.2.0.tar.gz"
)
local <- file.path(
    "annotations",
    "refdata-cellranger-mm10-1.2.0.tar.gz"
)
download.file(remote, local)
untar(local)
rm(remote, local)

bcb <- loadCellRanger(
    uploadDir = file.path("data", "cellranger"),
    refDataDir = file.path("annotations", "refdata-cellranger-mm10-1.2.0"),
    sampleMetadataFile = file.path("meta", "sample_metadata.xlsx"),
    interestingGroups = c("genotype", "age")
)
saveData(bcb, dir = file.path("data", Sys.Date()))
