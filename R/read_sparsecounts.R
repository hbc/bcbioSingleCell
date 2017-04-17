#' Read a MatrixMart file, setting the row and column names
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param matfile MatrixMart file to read
#'
#' @return a sparse matrix of counts
#' @export
read_sparsecounts <- function(matfile) {
    if (!file.exists(matfile)) {
        stop("count matrix could not be found")
    }
    if (file_ext(matfile) == "gz") {
        stem <- file_path_sans_ext(matfile)
    }
    else {
        stem <- matfile
    }
    rowfile <- paste0(stem, ".rownames")
    colfile <- paste0(stem, ".colnames")
    if (!file.exists(rowfile)) {
        stop("rownames file could not be found")
    }
    if (!file.exists(colfile)) {
        stop("rownames file could not be found")
    }
    counts <- readMM(matfile)
    rownames(counts) <- read_lines(rowfile)
    colnames(counts) <- read_lines(colfile)
    return(as(counts, "dgCMatrix"))
}
