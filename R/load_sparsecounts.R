#' Load a gzipped MatrixMart file, setting the row and column names
#'
#' @author Rory Kirchner
#'
#' @param gzfile gzipped MatrixMart file to load
#'
#' @return a sparse matrix of counts
#' @export
load_sparsecounts <- function(gzfile) {
    stem <- file_path_sans_ext(gzfile)
    rowfile <- paste0(stem, ".rownames")
    colfile <- paste0(stem, ".colnames")
    counts <- readMM(gzfile)
    rownames(counts) <- read_lines(rowfile)
    colnames(counts) <- read_lines(colfile)
    return(counts)
}
