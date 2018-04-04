# cellranger dataset
# 2018-03-25
library(devtools)
library(Matrix)
library(readr)
load_all()

# Include the top 500 genes (rows) and cells (columns)
uploadDir <- "inst/extdata/cellranger"
sampleDir <- file.path(
    uploadDir,
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
countsPerGene <- rowSums(mat) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(countsPerGene))
countsPerCell <- colSums(mat) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- sort(names(countsPerCell))
mat <- mat[genes, cells]

# Update the `genes.tsv` to match
match <- match(x = rownames(mat), table = gene2symbol[[1L]])
stopifnot(!any(is.na(match)))
gene2symbol <- gene2symbol[match, ]
stopifnot(identical(rownames(mat), gene2symbol[[1L]]))

# Write update files to disk
writeMM(mat, file = matFile)
write_tsv(gene2symbol, path = genesFile, col_names = FALSE)
write_lines(colnames(mat), path = barcodesFile)

# cellranger_small =============================================================
cellranger_small <- loadCellRanger(
    uploadDir = uploadDir,
    refdataDir = file.path(uploadDir, "refdata-cellranger-hg19-1.2.0"),
    sampleMetadataFile = file.path(uploadDir, "metadata.csv")
)

# Apply example filtering without excluding any cells
cellranger_small <- filterCells(
    cellranger_small,
    minUMIs = 0,
    minGenes = 0,
    maxMitoRatio = Inf,
    minNovelty = 0
)

use_data(cellranger_small, compress = "xz", overwrite = TRUE)
