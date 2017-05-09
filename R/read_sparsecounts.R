#' Read a MatrixMart file, setting the row and column names
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param matrix_file MatrixMart file to read
#'
#' @return a sparse matrix of counts
#' @export
read_counts <- function(matrix_file) {
    row_file <- paste0(matrix_file, ".rownames")
    col_file <- paste0(matrix_file, ".colnames")

    if (!file.exists(matrix_file)) {
        stop("MatrixMart counts file missing")
    }
    if (!file.exists(row_file)) {
        stop("Rownames file could not be found")
    }
    if (!file.exists(col_file)) {
        stop("Colnames file could not be found")
    }

    counts <- readMM(matrix_file)
    rownames(counts) <- read_lines(row_file)
    colnames(counts) <- read_lines(col_file)

    # Coerce dgTMatrix to dgCMatrix
    return(as(counts, "dgCMatrix"))
}
