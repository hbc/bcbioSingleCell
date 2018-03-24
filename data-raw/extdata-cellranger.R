# cellranger dataset
# 2018-03-23
library(devtools)
library(Matrix)
library(readr)
load_all()

# Include the top 500 genes (rows) and cells (columns)
cellrangerDir <- "inst/extdata/cellranger"
sampleDir <- file.path(
    cellrangerDir,
    "aggregation",
    "outs",
    "filtered_gene_bc_matrices",
    "hg19"
)

matFile <- file.path(sampleDir, "matrix.mtx")
genesFile <- file.path(sampleDir, "genes.tsv")
barcodesFile <- file.path(sampleDir, "barcodes.tsv")

# Prepare the sparse matrix
mat <- readMM(matFile)
gene2symbol <- read_tsv(
    file = genesFile,
    col_names = c("geneID", "geneName")
)
rownames <- gene2symbol[["geneID"]]
colnames <- read_lines(barcodesFile)
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

# Update the `genes.tsv` to match
match <- match(x = rownames(mat), table = gene2symbol[[1L]])
stopifnot(!any(is.na(match)))
barcodes <- barcodes[match, ]
stopifnot(identical(colnames(mat), barcodes[[1L]]))

# Write update files to disk
writeMM(mat, file = matFile)
write_lines(gene2symbol, path = genesFile)
write_lines(colnames(mat), path = barcodesFile)

# cellranger_small =============================================================
cellranger_small <- loadSingleCell(
    uploadDir = cellrangerDir,
    sampleMetadataFile = file.path(cellrangerDir, "metadata.csv"),
    organism = "Homo sapiens"
)

# Apply example filtering without excluding any cells
cellranger_small <- filterCells(
    cellranger_small,
    minUMIs = 0,
    minGenes = 0,
    maxMitoRatio = Inf,
    minNovelty = 0
)

saveData(cellranger_small, dir = "inst/extdata", compress = "xz")
