# nolint start
#
# Load bcbio Single-Cell RNA-Seq Data
#
# Michael Steinbaugh
# 2018-02-26
#
# Example script for *Mus musculus* with Ensembl 90 annotations
#
# Latest version of this script is available here:
# script <- system.file(
#     file.path("R_scripts", "load.R"),
#     package = "bcbioSingleCell")
# file.edit(script)
#
# nolint end

library(bcbioSingleCell)
bcb <- loadSingleCell(
    uploadDir = path("data", "indrop_rnaseq"),
    interestingGroups = c("genotype", "treatment"),
    sampleMetadataFile = path("meta", "sample_metadata.xlsx"),
    ensemblVersion = 90L
)
saveData(bcb, dir = path("data", Sys.Date()))
