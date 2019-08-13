## inDrops example data
## Using harvard-indrop-v3 barcodes.
## Updated 2019-08-12.

library(usethis)
library(pryr)
library(tidyverse)
library(Matrix)

## Restrict to 1 MB.
## Use `pryr::object_size` instead of `utils::object.size`.
limit <- structure(1e6, class = "object_size")

## Minimal example bcbio upload directory =======================================
## Include the top 500 genes (rows) and cells (columns).
upload_dir <- "inst/extdata/indrops"
sample <- "multiplexed-AAAAAAAA"

counts_file <- file.path(
    upload_dir,
    sample,
    paste0(sample, ".mtx")
)
rownames_file <- file.path(
    upload_dir,
    sample,
    paste0(sample, ".mtx.rownames")
)
colnames_file <- file.path(
    upload_dir,
    sample,
    paste0(sample, ".mtx.colnames")
)
barcodes_file <- file.path(
    upload_dir,
    sample,
    paste0(sample, "-barcodes.tsv")
)

stopifnot(all(file.exists(
    c(counts_file, rownames_file, colnames_file, barcodes_file)
)))

## Prepare the sparse matrix.
counts <- Matrix::readMM(counts_file)
rownames <- readr::read_lines(rownames_file)
colnames <- readr::read_lines(colnames_file)
stopifnot(
    identical(nrow(counts), length(rownames)),
    identical(ncol(counts), length(colnames))
)
rownames(counts) <- rownames
colnames(counts) <- colnames

## Subset the matrix to include only the top genes and cells.
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 50L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 100L)
cells <- sort(names(top_cells))

counts <- counts[genes, cells]

## Update the `barcodes.tsv` file to match.
barcodes <- readr::read_tsv(barcodes_file, col_names = FALSE)
match <- match(x = colnames(counts), table = barcodes[[1L]])
stopifnot(!any(is.na(match)))
barcodes <- barcodes[match, ]
stopifnot(identical(colnames(counts), barcodes[[1L]]))

## Write update files to disk.
Matrix::writeMM(counts, file = counts_file)
readr::write_lines(rownames(counts), path = rownames_file)
readr::write_lines(colnames(counts), path = colnames_file)
readr::write_tsv(barcodes, path = barcodes_file, col_names = FALSE)

## bcbioSingleCell object =======================================================
bcb <- bcbioSingleCell(
    uploadDir = upload_dir,
    sampleMetadataFile = file.path(upload_dir, "metadata.csv"),
    organism = "Homo sapiens",
    ensemblRelease = 90L
)

## Report the size of each slot in bytes.
lapply(coerceS4ToList(bcb), object_size)
object_size(bcb)
stopifnot(object_size(bcb) < limit)
validObject(bcb)

use_data(bcb, compress = "xz", overwrite = TRUE)
