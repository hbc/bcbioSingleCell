#' Read a gzipped MatrixMart file, setting the row and column names
#'
#' @author Rory Kirchner
#'
#' @keywords internal
#'
#' @param gzfile gzipped MatrixMart file to read
#'
#' @return a sparse matrix of counts
#' @export
read_sparsecounts_gz <- function(gzfile) {
    stem <- file_path_sans_ext(gzfile)
    rowfile <- paste0(stem, ".rownames")
    colfile <- paste0(stem, ".colnames")
    counts <- readMM(gzfile)
    rownames(counts) <- read_lines(rowfile)
    colnames(counts) <- read_lines(colfile)
    return(counts)
}
