# Cell Ranger example data
# 2018-11-14
# 4k PBMCs from a Healthy Donor
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k

library(pryr)
library(tidyverse)
library(Matrix)

# Restrict to 1 MB.
# Use `pryr::object_size()` instead of `utils::object.size()`.
limit <- structure(1e6, class = "object_size")

# Complete dataset =============================================================
upload_dir <- initDir("data-raw/cellranger")
# Example dataset contains a single sample ("pbmc4k").
outs_dir <- initDir(file.path(upload_dir, "pbmc", "outs"))

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
gtf_file <- file.path(upload_dir, "genes.gtf")
if (!file.exists(gtf_file)) {
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
        destfile = gtf_file
    )
}

# Note that this blows up in memory too much to run on RStudio AMI.
sce <- CellRanger(
    uploadDir = upload_dir,
    organism = "Homo sapiens",
    gffFile = gtf_file
)
object_size(sce)
assignAndSaveData(name = "pbmc", object = sce, dir = "data-raw")

# cellranger =============================================================
counts <- counts(sce)

# Subset the matrix to include only the top genes and cells.
top_genes <- Matrix::rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 50L)
genes <- sort(names(top_genes))

top_cells <- Matrix::colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 100L)
cells <- sort(names(top_cells))

# Subset the original pbmc dataset to contain only top genes and cells.
sce <- sce[genes, cells]

# Include only minimal metadata columns in rowRanges.
mcols <- mcols(rowRanges(sce))
mcols <- mcols[, c("broadClass", "geneBiotype", "geneID", "geneName")]
# TODO Consider making this Rle factor level step a function in basejump.
mcols <- DataFrame(lapply(
    X = mcols,
    FUN = function(x) {
        if (is(x, "Rle")) {
            x <- decode(x)
            if (is.factor(x)) {
                x <- droplevels(x)
            }
            Rle(x)
        } else {
            I(x)
        }
    }
))
mcols(rowRanges(sce)) <- mcols

# Report the size of each slot in bytes.
vapply(
    X = coerceS4ToList(sce),
    FUN = object_size,
    FUN.VALUE = numeric(1L)
)
object_size(sce)
stopifnot(object_size(sce) < limit)
stopifnot(validObject(sce))

cellranger <- sce
usethis::use_data(cellranger, compress = "xz", overwrite = TRUE)



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
counts <- Matrix::readMM(file = file.path(input_dir, "matrix.mtx"))
counts <- counts[seq_len(100), seq_len(100)]
Matrix::writeMM(counts, file = file.path(output_dir, "matrix.mtx"))

genes <- readr::read_tsv(
    file = file.path(input_dir, "genes.tsv"),
    col_names = c("geneID", "geneName"),
    n_max = 100
)
readr::write_tsv(
    x = genes,
    path = file.path(output_dir, "genes.tsv"),
    col_names = FALSE
)

barcodes <- readr::read_lines(
    file = file.path(input_dir, "barcodes.tsv"),
    n_max = 100
)
readr::write_lines(
    x = barcodes,
    path = file.path(output_dir, "barcodes.tsv")
)
