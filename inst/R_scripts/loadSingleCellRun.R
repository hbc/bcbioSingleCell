# Example script for Mus musculus with Ensembl 88 annotations
library(bcbioSingleCell)
dir.create("annotations", showWarnings = FALSE)
download.file(
    file.path("http://ftp.ensembl.org",
              "pub",
              "release-88",
              "mus_musculus",
              "Mus_musculus.GRCm38.88.chr_patch_hapl_scaff.gtf.gz"),
    file.path("annotations",
              "Mus_musculus.GRCm38.88.chr_patch_hapl_scaff.gtf.gz"))
bcb <- loadSingleCellRun(
    file.path("data", "indrop_rnaseq"),
    sampleMetadataFile = file.path("meta", "sample_metadata.xlsx"),
    interestingGroups = c("genotype", "treatment"),
    ensemblVersion = 88,
    gtfFile = file.path("annotations",
                        "Mus_musculus.GRCm38.88.chr_patch_hapl_scaff.gtf.gz"))
saveData(bcb)
