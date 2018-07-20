# indrops_small (harvard-indrop-v3)
# 2018-07-12

library(assertive)
library(devtools)
library(tidyverse)
library(Matrix)
load_all()

# Include the top 500 genes (rows) and cells (columns)
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
assert_all_are_existing_files(
    c(counts_file, rownames_file, colnames_file, barcodes_file)
)

# Prepare the sparse matrix
counts <- readMM(counts_file)
rownames <- read_lines(rownames_file)
colnames <- read_lines(colnames_file)
stopifnot(identical(nrow(counts), length(rownames)))
stopifnot(identical(ncol(counts), length(colnames)))
rownames(counts) <- rownames
colnames(counts) <- colnames

# Subset the matrix to include only the top genes and cells
top_genes <- rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- sort(names(top_cells))

counts <- counts[genes, cells]

# Update the `barcodes.tsv` file to match
barcodes <- read_tsv(barcodes_file, col_names = FALSE)
match <- match(x = colnames(counts), table = barcodes[[1L]])
stopifnot(!any(is.na(match)))
barcodes <- barcodes[match, ]
stopifnot(identical(colnames(counts), barcodes[[1L]]))

# Write update files to disk
writeMM(counts, file = counts_file)
write_lines(rownames(counts), path = rownames_file)
write_lines(colnames(counts), path = colnames_file)
write_tsv(barcodes, path = barcodes_file, col_names = FALSE)

# indrops_small ================================================================
indrops_small <- bcbioSingleCell(
    uploadDir = upload_dir,
    sampleMetadataFile = file.path(upload_dir, "metadata.csv"),
    organism = "Homo sapiens",
    ensemblRelease = 90L
)

# Apply example filtering without excluding any cells
indrops_small <- filterCells(
    object = indrops_small,
    minUMIs = 0,
    minGenes = 0,
    minNovelty = 0,
    maxMitoRatio = 1,
    minCellsPerGene = 0
)

# Require 500 cells, 500 genes
assert_are_identical(dim(indrops_small), c(500L, 500L))

use_data(indrops_small, compress = "xz", overwrite = TRUE)
