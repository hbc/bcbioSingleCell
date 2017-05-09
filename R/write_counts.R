#' Write MatrixMart .mtx, .colnames, and .rownames files
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param counts Sparse counts matrix
#' @param gzip Compress matrix file with gzip
#'
#' @return Path to matrix file
#' @export
write_counts <- function(counts, gzip = FALSE) {
    # Write MatrixMart .mtx file
    matrix_file <- paste0(deparse(substitute(counts)), ".mtx")
    writeMM(counts, matrix_file)

    # Write barcodes (colnames)
    barcodes <- colnames(counts)
    barcodes_file <- paste0(matrix_file, ".colnames")
    write_lines(barcodes, barcodes_file)

    # Write gene names (rownames)
    genes <- rownames(counts)
    genes_file <- paste0(matrix_file, ".rownames")
    write_lines(genes, genes_file)

    # gzip the matrix, if desired
    if (isTRUE(gzip)) {
        matrix_file <- gzip(matrix_file)
    }

    return(matrix_file)
}
