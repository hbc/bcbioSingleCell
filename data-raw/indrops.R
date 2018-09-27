# inDrops example data
# Using harvard-indrop-v3 barcodes
# 2018-09-27

library(tidyverse)
library(Seurat)
library(Matrix)

# Restrict to 1 MB per file.
limit <- structure(1e6, class = "object_size")

# Minimal example bcbio upload directory =======================================
# Include the top 500 genes (rows) and cells (columns).
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

# Prepare the sparse matrix.
counts <- readMM(counts_file)
rownames <- read_lines(rownames_file)
colnames <- read_lines(colnames_file)
stopifnot(identical(nrow(counts), length(rownames)))
stopifnot(identical(ncol(counts), length(colnames)))
rownames(counts) <- rownames
colnames(counts) <- colnames

# Subset the matrix to include only the top genes and cells.
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- sort(names(top_cells))

counts <- counts[genes, cells]

# Update the `barcodes.tsv` file to match.
barcodes <- read_tsv(barcodes_file, col_names = FALSE)
match <- match(x = colnames(counts), table = barcodes[[1L]])
stopifnot(!any(is.na(match)))
barcodes <- barcodes[match, ]
stopifnot(identical(colnames(counts), barcodes[[1L]]))

# Write update files to disk.
writeMM(counts, file = counts_file)
write_lines(rownames(counts), path = rownames_file)
write_lines(colnames(counts), path = colnames_file)
write_tsv(barcodes, path = barcodes_file, col_names = FALSE)

# bcbioSingleCell object =======================================================
sce <- bcbioSingleCell(
    uploadDir = upload_dir,
    sampleMetadataFile = file.path(upload_dir, "metadata.csv"),
    organism = "Homo sapiens",
    ensemblRelease = 90L
)

# Include only minimal metadata columns in rowRanges.
mcols(rowRanges(sce)) <- mcols(rowRanges(sce)) %>%
    as("tbl_df") %>%
    select(rowname, broadClass, geneBiotype, geneID, geneName) %>%
    mutate_if(is.factor, droplevels) %>%
    as("DataFrame")

# Report the size of each slot in bytes.
vapply(
    X = coerceS4ToList(sce),
    FUN = object.size,
    FUN.VALUE = numeric(1L)
)
stopifnot(object.size(sce) < limit)
stopifnot(validObject(sce))

indrops_small <- sce
devtools::use_data(indrops_small, compress = "xz", overwrite = TRUE)
