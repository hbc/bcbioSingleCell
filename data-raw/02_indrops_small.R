# inDrops example dataset (indrops_small)
# Using harvard-indrop-v3 barcodes
# 2018-07-26

library(devtools)
library(tidyverse)
library(Seurat)
library(Matrix)
load_all()

# Minimal example bcbio upload directory =======================================
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

stopifnot(all(file.exists(
    c(counts_file, rownames_file, colnames_file, barcodes_file)
)))

# Prepare the sparse matrix
counts <- readMM(counts_file)
rownames <- read_lines(rownames_file)
colnames <- read_lines(colnames_file)
stopifnot(identical(nrow(counts), length(rownames)))
stopifnot(identical(ncol(counts), length(colnames)))
rownames(counts) <- rownames
colnames(counts) <- colnames

# Subset the matrix to include only the top genes and cells
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
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

# bcbioRNASeq object ===========================================================
bcb <- bcbioSingleCell(
    uploadDir = upload_dir,
    sampleMetadataFile = file.path(upload_dir, "metadata.csv"),
    organism = "Homo sapiens",
    ensemblRelease = 90L
)

# Apply example filtering without excluding any cells
bcb <- filterCells(
    object = bcb,
    minUMIs = 0,
    maxUMIs = Inf,
    minGenes = 0,
    maxGenes = Inf,
    minNovelty = 0,
    maxMitoRatio = 1,
    minCellsPerGene = 0
)

# Require 500 cells, 500 genes
stopifnot(identical(dim(indrops_small), c(500L, 500L)))

# seurat =======================================================================
# Let's handoff to seurat to perform dimensionality reduction and clustering,
# then slot the DR data in our bcbioRNASeq object
seurat <- as(bcb, "seurat") %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(
        resolution = seq(from = 0.4, to = 1.2, by = 0.4)
    ) %>%
    RunTSNE() %>%
    # Requires `umap-learn` Python package.
    # Install with conda or pip.
    RunUMAP() %>%
    SetAllIdent(id = "res.0.4")

# Coerce seurat to SingleCellExperiment, which contains `reducedDims` slot
seurat_sce <- as(seurat, "SingleCellExperiment")
stopifnot(identical(
    names(reducedDims(seurat_sce)),
    c("PCA", "TSNE", "UMAP")
))
stopifnot(identical(dimnames(bcb), dimnames(seurat_sce)))



# Save =========================================================================
# Slot the reduced dimensions into our bcbioSingleCell object
reducedDims(bcb) <- reducedDims(seurat_sce)
indrops_small <- bcb
use_data(indrops_small, compress = "xz", overwrite = TRUE)
