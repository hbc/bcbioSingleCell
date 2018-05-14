# indrops_small (harvard-indrop-v3)
# 2018-05-14
library(devtools)
load_all()
library(tidyverse)
library(Matrix)

# Include the top 500 genes (rows) and cells (columns)
uploadDir <- "inst/extdata/indrops"
sampleName <- "multiplexed-AAAAAAAA"

matFile <- file.path(
    uploadDir,
    sampleName,
    paste0(sampleName, ".mtx")
)
rownamesFile <- file.path(
    uploadDir,
    sampleName,
    paste0(sampleName, ".mtx.rownames")
)
colnamesFile <- file.path(
    uploadDir,
    sampleName,
    paste0(sampleName, ".mtx.colnames")
)
barcodesFile <- file.path(
    uploadDir,
    sampleName,
    paste0(sampleName, "-barcodes.tsv")
)

# Prepare the sparse matrix
mat <- readMM(matFile)
rownames <- read_lines(rownamesFile)
colnames <- read_lines(colnamesFile)
stopifnot(identical(nrow(mat), length(rownames)))
stopifnot(identical(ncol(mat), length(colnames)))
rownames(mat) <- rownames
colnames(mat) <- colnames

# Subset the matrix to include only the top genes and cells
countsPerGene <- rowSums(mat) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(countsPerGene))
countsPerCell <- colSums(mat) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- sort(names(countsPerCell))
mat <- mat[genes, cells]

# Update the `barcodes.tsv` file to match
barcodes <- read_tsv(barcodesFile, col_names = FALSE)
match <- match(x = colnames(mat), table = barcodes[[1L]])
stopifnot(!any(is.na(match)))
barcodes <- barcodes[match, ]
stopifnot(identical(colnames(mat), barcodes[[1L]]))

# Write update files to disk
writeMM(mat, file = matFile)
write_lines(rownames(mat), path = rownamesFile)
write_lines(colnames(mat), path = colnamesFile)
write_tsv(barcodes, path = barcodesFile, col_names = FALSE)

# indrops_small ================================================================
indrops_small <- bcbioSingleCell(
    uploadDir = uploadDir,
    sampleMetadataFile = file.path(uploadDir, "metadata.csv"),
    organism = "Homo sapiens",
    ensemblRelease = 90L
)

# Apply example filtering without excluding any cells
indrops_small <- filterCells(
    object = indrops_small,
    minUMIs = 0,
    minGenes = 0,
    minNovelty = 0,
    maxMitoRatio = Inf,
    minCellsPerGene = 0
)

use_data(indrops_small, compress = "xz", overwrite = TRUE)
