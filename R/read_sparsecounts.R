#' Read a MatrixMart file, setting the row and column names
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param mat_file MatrixMart file to read
#'
#' @return a sparse matrix of counts
#' @export
read_sparsecounts <- function(mat_file) {
    if (!file.exists(mat_file)) {
        stop("count matrix could not be found")
    }
    if (file_ext(mat_file) == "gz") {
        stem <- file_path_sans_ext(mat_file)
    }
    else {
        stem <- mat_file
    }
    row_file <- paste0(stem, ".rownames")
    col_file <- paste0(stem, ".colnames")
    if (!file.exists(row_file)) {
        stop("rownames file could not be found")
    }
    if (!file.exists(col_file)) {
        stop("rownames file could not be found")
    }
    counts <- readMM(mat_file)
    rownames(counts) <- read_lines(row_file)
    colnames(counts) <- read_lines(col_file)
    return(as(counts, "dgCMatrix"))
}
