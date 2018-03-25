# harvard_indrop_v3 dataset
# 2018-03-25
library(devtools)
library(Matrix)
library(readr)
load_all()

# Include the top 500 genes (rows) and cells (columns)
indropDir <- "inst/extdata/harvard_indrop_v3"
sampleName <- "multiplexed-AAAAAAAA"

matFile <- file.path(
    indropDir,
    sampleName,
    paste0(sampleName, ".mtx")
)
rownamesFile <- file.path(
    indropDir,
    sampleName,
    paste0(sampleName, ".mtx.rownames")
)
colnamesFile <- file.path(
    indropDir,
    sampleName,
    paste0(sampleName, ".mtx.colnames")
)
barcodesFile <- file.path(
    indropDir,
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
countsPerGene <- Matrix::rowSums(mat) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- names(countsPerGene)
countsPerCell <- Matrix::colSums(mat) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- names(countsPerCell)
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

# indrop_small =================================================================
indrop_small <- loadSingleCell(
    uploadDir = indropDir,
    sampleMetadataFile = file.path(indropDir, "metadata.csv"),
    organism = "Homo sapiens"
)

# Apply example filtering without excluding any cells
indrop_small <- filterCells(
    indrop_small,
    minUMIs = 0,
    minGenes = 0,
    maxMitoRatio = Inf,
    minNovelty = 0
)

saveData(indrop_small, dir = "inst/extdata", compress = "xz")
