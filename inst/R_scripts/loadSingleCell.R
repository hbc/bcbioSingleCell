# Latest version of this script is available here:
# script <- system.file(
#     file.path("R_scripts", "loadSingleCell.R"),
#     package = "bcbioSingleCell")
# file.edit(script)

# Example script for Mus musculus with Ensembl 88 annotations
library(bcbioSingleCell)

dir.create("annotations", showWarnings = FALSE)
download.file(
    url = file.path(
        "http://ftp.ensembl.org",
        "pub",
        "release-88",
        "mus_musculus",
        "Mus_musculus.GRCm38.88.chr_patch_hapl_scaff.gtf.gz"),
    destfile = file.path(
        "annotations",
        "Mus_musculus.GRCm38.88.chr_patch_hapl_scaff.gtf.gz"))

bcb <- loadSingleCell(
    uploadDir = file.path("data", "indrop_rnaseq"),
    interestingGroups = c("genotype", "treatment"),
    sampleMetadataFile = file.path("meta", "sample_metadata.xlsx"),
    gtfFile = file.path(
        "annotations",
        "Mus_musculus.GRCm38.88.chr_patch_hapl_scaff.gtf.gz"),
    ensemblVersion = 88)

# Back up all data inside bcbio object
flatFiles <- flatFiles(bcb)

saveData(bcb, flatFiles, dir = "data")
