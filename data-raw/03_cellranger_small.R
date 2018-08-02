# Cell Ranger Example Data
# 2018-08-01
# 4k PBMCs from a Healthy Donor
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k

library(devtools)
library(tidyverse)
library(Seurat)
library(Matrix)
load_all()

# Complete dataset =============================================================
# Directory structure:
# - pbmc
# - outs
# - filtered_gene_bc_matrices
# - GRCh38
upload_dir <- "data-raw/cellranger"
unlink(upload_dir, recursive = TRUE)
dir.create(upload_dir, recursive = TRUE)
# Example dataset contains a single sample ("pbmc4k")
outs_dir <- file.path(upload_dir, "pbmc", "outs")
dir.create(outs_dir, recursive = TRUE, showWarnings = FALSE)
tar_file <- file.path(upload_dir, "pbmc.tar.gz")
if (!file.exists(tar_file)) {
    download.file(
        url = paste(
            "http://cf.10xgenomics.com",
            "samples",
            "cell-exp",
            "2.1.0",
            "pbmc4k",
            "pbmc4k_filtered_gene_bc_matrices.tar.gz",
            sep = "/"
        ),
        destfile = tar_file
    )
}
untar(tarfile = tar_file, exdir = outs_dir)
stopifnot(identical(dir(outs_dir), "filtered_gene_bc_matrices"))

# Ensembl 92, expect 483 unannotated genes as warning
# dim(sce)
# [1] 33694  4340
# Too large to save inside package
sce <- suppressWarnings(
    readCellRanger(uploadDir = upload_dir, organism = "Homo sapiens")
)
saveData(sce, dir = "~")



# Minimal SingleCellExperiment =================================================
# Subset to only include the top 500 genes (rows) and cells (columns)
counts <- counts(sce)

# Subset the matrix to include only the top genes and cells
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- sort(names(top_cells))

# Subset the full PBMC dataset to be minimal
sce <- sce[genes, cells]



# Minimal example cellranger directory =========================================
input_dir <- file.path(
    upload_dir,
    "pbmc",
    "outs",
    "filtered_gene_bc_matrices",
    "GRCh38"
)
stopifnot(dir.exists(input_dir))
output_dir <- file.path(
    "inst",
    "extdata",
    "cellranger",
    "pbmc",
    "outs",
    "filtered_gene_bc_matrices",
    "GRCh38"
)
unlink("inst/extdata/cellranger", recursive = TRUE)
dir.create(output_dir, recursive = TRUE)

# Prepare the sparse matrix
counts <- readMM(file.path(input_dir, "matrix.mtx"))
counts <- counts[1:500, 1:500]
writeMM(counts, file = file.path(output_dir, "matrix.mtx"))

genes <- read_tsv(
    file = file.path(input_dir, "genes.tsv"),
    col_names = c("geneID", "geneName"),
    n_max = 500
)
write_tsv(
    x = genes,
    path = file.path(output_dir, "genes.tsv"),
    col_names = FALSE
)

barcodes <- read_lines(
    file = file.path(input_dir, "barcodes.tsv"),
    n_max = 500
)
write_lines(
    x = barcodes,
    path = file.path(output_dir, "barcodes.tsv")
)



# seurat_small =================================================================
# Let's handoff to seurat to perform dimensionality reduction and clustering,
# then slot the DR data in our bcbioRNASeq object
seurat_small <- as(sce, "seurat") %>%
    convertGenesToSymbols() %>%
    NormalizeData() %>%
    FindVariableGenes(do.plot = FALSE) %>%
    ScaleData() %>%
    RunPCA(do.print = FALSE) %>%
    FindClusters(resolution = seq(from = 0.4, to = 1.2, by = 0.4)) %>%
    RunTSNE() %>%
    # Requires Python `umap-learn` package
    RunUMAP() %>%
    SetAllIdent(id = "res.0.4")



# all_markers_small ============================================================
all_markers_small <- seurat_small %>%
    FindAllMarkers() %>%
    # Sanitize, including more robust gene annotation information
    sanitizeSeuratMarkers(rowRanges = rowRanges(seurat_small))



# known_markers_small ==========================================================
known_markers_small <- knownMarkersDetected(
    all = all_markers_small,
    known = cell_type_markers[["homoSapiens"]]
)



# cellranger_small =============================================================
# Convert rows (geneName) back to Ensembl IDs (geneID)
seurat_sce <- seurat_small %>%
    as("SingleCellExperiment") %>%
    convertSymbolsToGenes()
# Ensure that dimensional reduction data is slotted correctly
stopifnot(identical(
    names(reducedDims(seurat_sce)),
    c("PCA", "TSNE", "UMAP")
))
stopifnot(identical(dimnames(sce), dimnames(seurat_sce)))
colData(sce) <- colData(seurat_sce)
reducedDims(sce) <- reducedDims(seurat_sce)
cellranger_small <- sce



# Save =========================================================================
use_data(
    cellranger_small,
    seurat_small,
    all_markers_small,
    known_markers_small,
    compress = "xz",
    overwrite = TRUE
)
