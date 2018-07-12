# Cell Ranger Example Data
# 2018-07-12
# 4k PBMCs from a Healthy Donor
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k

library(assertive)
library(devtools)
library(tidyverse)
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
untar(
    tarfile = tar_file,
    exdir = outs_dir
)
stopifnot(identical(dir(outs_dir), "filtered_gene_bc_matrices"))

pbmc <- readCellRanger(
    uploadDir = upload_dir,
    organism = "Homo sapiens"
)
# Ensembl 92, expect 483 unannotated genes
# dim(pbmc)
# [1] 33694  4340
# Too large to save inside package
saveData(pbmc, dir = "~")



# cellranger_small =============================================================
# Subset to only include the top 500 genes (rows) and cells (columns)
counts <- counts(pbmc)

# Subset the matrix to include only the top genes and cells
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- sort(names(top_cells))

cellranger_small <- pbmc[genes, cells]
use_data(cellranger_small, compress = "xz", overwrite = TRUE)



# readCellRanger extdata example ===============================================
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
