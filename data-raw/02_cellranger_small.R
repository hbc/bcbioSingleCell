# cellranger_small
# 2018-06-04
# 4k PBMCs from a Healthy Donor
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k
library(devtools)
library(tidyverse)
library(Matrix)
load_all()

unlink("data-raw/pbmc4k", recursive = TRUE)
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
    destfile = "data-raw/pbmc4k.tar.gz"
)
untar(
    tarfile = "data-raw/pbmc4k.tar.gz",
    exdir = "data-raw/pbmc4k"
)

pbmc4k <- readCellRanger("data-raw/pbmc4k")
# Ensembl 92, 483 unannotated genes
# dim(pbmc4k)
# [1] 33694  4340
# Too large to save inside package
saveData(pbmc4k, dir = "~")

# Subset to only include the top 500 genes (rows) and cells (columns)
counts <- counts(pbmc4k)

# Subset the matrix to include only the top genes and cells
top_genes <- rowSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
genes <- sort(names(top_genes))

top_cells <- colSums(counts) %>%
    sort(decreasing = TRUE) %>%
    head(n = 500L)
cells <- sort(names(top_cells))

cellranger_small <- pbmc4k[genes, cells]

use_data(cellranger_small, compress = "xz", overwrite = TRUE)
