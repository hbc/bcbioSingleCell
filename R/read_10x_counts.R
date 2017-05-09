#' Read 10x Genomics counts data
#'
#' Generate a sparse matrix from barcodes.tsv, genes.tsv, and matrix.mtx files
#'
#' @author Michael Steinbaugh
#'
#' @param data_dir Data directory cont
#'
#' @return Sparse counts matrix (dgTMatrix)
#' @export
read_10x_counts <- function(data_dir) {
    barcodes_file <- file.path(data_dir, "barcodes.tsv")
    genes_file <- file.path(data_dir, "genes.tsv")
    matrix_file <- file.path(data_dir, "matrix.mtx")

    if (!file.exists(barcodes_file)) {
        stop("Barcodes file missing")
    }
    if (!file.exists(genes_file)) {
        stop("Gene names file missing")
    }
    if (!file.exists(matrix_file)) {
        stop("MatrixMart counts matrix file missing")
    }

    # Read MatrixMart expression counts matrix
    sparsecounts <- readMM(matrix_file)

    # Assign barcodes to colnames
    sparse_cols <- readLines(barcodes_file)
    # Check for `-1` suffix in barcode names, indicative of 10X data
    pattern <- "\\-1$"
    if (all(grepl(pattern, sparse_cols))) {
        sparse_cols <- gsub(pattern, "", sparse_cols)
    }
    colnames(sparsecounts) <- sparse_cols

    # Assign gene names (symbols) to rownames
    sparse_rows <- read_tsv(genes_file,
                            col_names = c("ensembl_gene_id",
                                          "external_gene_name"))
    rownames(sparsecounts) <- make.unique(sparse_rows$external_gene_name)

    return(sparsecounts)
}
