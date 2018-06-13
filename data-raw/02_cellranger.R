# Cell Ranger Example Data
# 2018-06-13
# 4k PBMCs from a Healthy Donor
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k

library(assertive)
library(devtools)
library(tidyverse)
library(Matrix)
load_all()



# Complete dataset =============================================================
tar_file <- "data-raw/pbmc4k.tar.gz"
upload_dir <- "data-raw/pbmc4k"
unlink(upload_dir, recursive = TRUE)
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
    exdir = upload_dir
)

pbmc4k <- readCellRanger(uploadDir = upload_dir)
# Ensembl 92, 483 unannotated genes
# dim(pbmc4k)
# [1] 33694  4340
# Too large to save inside package
saveData(pbmc4k, dir = "~")



# cellranger_small =============================================================
# Subset to only include the top 500 genes (rows) and cells (columns)
counts <- counts(pbmc4k)

# Subset the matrix to include only the top genes and cells
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- sort(names(top_cells))

cellranger_small <- pbmc4k[genes, cells]
use_data(cellranger_small, compress = "xz", overwrite = TRUE)



# readCellRanger extdata example ===============================================
input_dir <- file.path(
    upload_dir,
    "filtered_gene_bc_matrices",
    "GRCh38"
)
output_dir <- file.path(
    "inst",
    "extdata",
    "cellranger",
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
