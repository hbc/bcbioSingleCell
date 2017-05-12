#' Read \href{https://www.10xgenomics.com/software/}{10x Genomics Chromium} run
#' output.
#'
#' @rdname read_10x
#' @author Michael Steinbaugh
#'
#' @param data_dir Data directory



#' @rdname read_10x
#' @description Read into sparse matrix from barcodes.tsv, genes.tsv, and
#'   matrix.mtx files.
#' @return Sparse counts matrix (dgCMatrix).
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
        stop("MatrixMart counts file missing")
    }

    # Read MatrixMart expression counts matrix
    counts <- readMM(matrix_file)

    # Assign barcodes to colnames
    barcodes <- readLines(barcodes_file)
    # Check for `-1` suffix in barcode names, indicative of 10X data
    pattern <- "\\-1$"
    if (all(grepl(pattern, barcodes))) {
        barcodes <- gsub(pattern, "", barcodes)
    }
    colnames(counts) <- barcodes

    # Assign gene names (symbols) to rownames
    genes <- read_tsv(genes_file,
                      col_names = c("ensembl_gene_id",
                                    "external_gene_name"))
    rownames(counts) <- make.unique(genes$external_gene_name)

    # Coerce dgTMatrix to dgCMatrix
    return(as(counts, "dgCMatrix"))
}
