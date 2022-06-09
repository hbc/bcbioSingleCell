## inDrops example data
## Using harvard-indrop-v3 barcodes.
## Updated 2022-06-09.
## nolint start
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(pipette)
})
## nolint end
load_all()
limit <- structure(2e6L, class = "object_size")
## Minimal example bcbio upload directory.
## Include the top 500 genes (rows) and cells (columns).
uploadDir <- file.path("..", "inst", "extdata", "indrops")
sample <- "multiplexed-AAAAAAAA"
countsFile <- file.path(
    uploadDir,
    sample,
    paste0(sample, ".mtx")
)
rownamesFile <- file.path(
    uploadDir,
    sample,
    paste0(sample, ".mtx.rownames")
)
colnamesFile <- file.path(
    uploadDir,
    sample,
    paste0(sample, ".mtx.colnames")
)
barcodesFile <- file.path(
    uploadDir,
    sample,
    paste0(sample, "-barcodes.tsv")
)
stopifnot(all(file.exists(
    c(countsFile, rownamesFile, colnamesFile, barcodesFile)
)))
barcodes <- import(barcodesFile, colnames = FALSE)
export(object = barcodes, con = barcodesFile, colnames = FALSE)
counts <- import(countsFile)
topGenes <-
    counts |>
    Matrix::rowSums() |>
    sort(decreasing = TRUE) |>
    head(n = 50L)
genes <- sort(names(topGenes))
cells <- barcodes[[1L]]
counts <- counts[genes, cells]
export(object = counts, con = countsFile)
## Create bcbioSingleCell object.
bcb <- bcbioSingleCell(
    uploadDir = uploadDir,
    sampleMetadataFile = file.path(uploadDir, "metadata.csv"),
    organism = "Homo sapiens",
    ensemblRelease = 90L
)
stopifnot(
    object.size(bcb) < limit,
    validObject(bcb)
)
use_data(bcb, compress = "xz", overwrite = TRUE)
