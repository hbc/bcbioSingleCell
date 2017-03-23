##' Write a gzipped MatrixMart file, storing the row and column names.
##'
##' @param sparse sparse matrix of counts
##' @param gzfile outfile with gzip extension
##' @return gzfile
##' @importFrom Matrix writeMM
##' @importFrom tools file_path_sans_ext
##' @importFrom R.utils gzip
##' @importFrom readr write_lines
##' @author Rory Kirchner
##' @export
write_sparsecounts <- function(sparse, gzfile) {
    if (file.exists(gzfile)) {
        return(gzfile)
    }
    stem <- tools::file_path_sans_ext(gzfile)
    rows <- rownames(sparse)
    cols <- colnames(sparse)
    rowfile <- paste0(stem, ".rownames")
    colfile <- paste0(stem, ".colnames")
    Matrix::writeMM(sparse, stem)
    R.utils::gzip(stem)
    readr::write_lines(rows, rowfile)
    readr::write_lines(cols, colfile)
    return(gzfile)
}
