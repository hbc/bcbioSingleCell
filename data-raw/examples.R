library(devtools)
load_all()

extdataDir <- system.file("extdata", package = "bcbioSingleCell")
uploadDir <- file.path(extdataDir, "harvard_indrop_v3")
sampleMetadataFile <- file.path(extdataDir, "harvard_indrop_v3.xlsx")
annotable <- annotable("Homo sapiens", release = 90)
bcb <- loadSingleCell(
    uploadDir = uploadDir,
    sampleMetadataFile = sampleMetadataFile,
    annotable = annotable)

examples <- list(
    bcb = bcb
)
use_data(examples, compress = "xz", overwrite = TRUE)
