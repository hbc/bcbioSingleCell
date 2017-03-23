##' Load a gzipped MatrixMart file, setting the row and column names.
##'
##' @param gzfile gzipped MatrixMart file to load
##' @return a sparse matrix of counts
##' @importFrom Matrix readMM
##' @importFrom tools file_path_sans_ext
##' @importFrom readr read_lines
##' @author Rory Kirchner
##' @export
load_sparsecounts <- function(gzfile) {
    stem <- tools::file_path_sans_ext(gzfile)
    rowfile <- paste0(stem, ".rownames")
    colfile <- paste0(stem, ".colnames")
    counts <- Matrix::readMM(gzfile)
    rownames(counts) <- readr::read_lines(rowfile)
    colnames(counts) <- readr::read_lines(colfile)
    return(counts)
}
