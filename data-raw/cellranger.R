# Cell Ranger example data
# 2018-09-19
# 4k PBMCs from a Healthy Donor
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k

library(tidyverse)
library(Seurat)
library(Matrix)

# Complete dataset =============================================================
upload_dir <- initializeDirectory("data-raw/cellranger")
# Example dataset contains a single sample ("pbmc4k").
outs_dir <- initializeDirectory(file.path(upload_dir, "pbmc", "outs"))

# Directory structure:
# - pbmc
# - outs
# - filtered_gene_bc_matrices
# - GRCh38
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

# Using Ensembl 84 GTF annotations.
gff_file <- file.path(upload_dir, "genes.gtf")
if (!file.exists(gff_file)) {
    download.file(
        url = paste(
            "ftp://ftp.ensembl.org",
            "pub",
            "release-84",
            "gtf",
            "homo_sapiens",
            "Homo_sapiens.GRCh38.84.gtf.gz",
            sep = "/"
        ),
        destfile = gff_file
    )
}

sce <- CellRanger(
    uploadDir = upload_dir,
    organism = "Homo sapiens",
    gffFile = gff_file
)

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

# Prepare the sparse matrix.
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

# cellranger_small =============================================================
counts <- counts(sce)

# Subset the matrix to include only the top genes and cells.
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- sort(names(top_cells))

cellranger_small <- sce[genes, cells]
devtools::use_data(cellranger_small, compress = "xz", overwrite = TRUE)
