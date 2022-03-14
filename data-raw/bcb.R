## inDrops example data
## Using harvard-indrop-v3 barcodes.
## Updated 2022-03-14.
suppressPackageStartupMessages({
    library(devtools)
    library(usethis)
    library(pryr)
    library(tidyverse)
    library(Matrix)
})
load_all()
## Restrict to 1 MB.
## Use `pryr::object_size` instead of `utils::object.size`.
limit <- structure(1e6L, class = "object_size")
## Minimal example bcbio upload directory.
## Include the top 500 genes (rows) and cells (columns).
uploadDir <- "inst/extdata/indrops"
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
## Prepare the sparse matrix.
counts <- Matrix::readMM(countsFile)
rownames <- readr::read_lines(rownamesFile)
colnames <- readr::read_lines(colnamesFile)
stopifnot(
    identical(nrow(counts), length(rownames)),
    identical(ncol(counts), length(colnames))
)
rownames(counts) <- rownames
colnames(counts) <- colnames
## Subset the matrix to include only the top genes and cells.
topGenes <-
    counts |>
    Matrix::rowSums() |>
    sort(decreasing = TRUE) |>
    head(n = 50L)
genes <- sort(names(topGenes))
topCells <-
    counts |>
    Matrix::colSums() |>
    sort(decreasing = TRUE) %>%
    head(n = 100L)
cells <- sort(names(topCells))
counts <- counts[genes, cells]
## Update the `barcodes.tsv` file to match.
barcodes <- readr::read_tsv(barcodesFile, col_names = FALSE)
match <- match(x = colnames(counts), table = barcodes[[1L]])
stopifnot(!any(is.na(match)))
barcodes <- barcodes[match, ]
stopifnot(identical(colnames(counts), barcodes[[1L]]))
## Write update files to disk.
Matrix::writeMM(counts, file = countsFile)
readr::write_lines(rownames(counts), path = rownamesFile)
readr::write_lines(colnames(counts), path = colnamesFile)
readr::write_tsv(barcodes, path = barcodesFile, col_names = FALSE)
## Create bcbioSingleCell object.
bcb <- bcbioSingleCell(
    uploadDir = uploadDir,
    sampleMetadataFile = file.path(uploadDir, "metadata.csv"),
    organism = "Homo sapiens",
    ensemblRelease = 90L
)
## Report the size of each slot in bytes.
lapply(coerceToList(bcb), object_size)
object_size(bcb)
stopifnot(object_size(bcb) < limit)
validObject(bcb)
use_data(bcb, compress = "xz", overwrite = TRUE)
