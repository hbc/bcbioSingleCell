#' Write a gzipped MatrixMart file, storing the row and column names.
#'
#' @author Rory Kirchner
#'
#' @param sparse sparse matrix of counts
#' @param gzfile outfile with gzip extension
#'
#' @return gzfile
#' @export
write_sparsecounts <- function(sparse, gzfile) {
    if (file.exists(gzfile)) {
        return(gzfile)
    }
    stem <- file_path_sans_ext(gzfile)
    rows <- rownames(sparse)
    cols <- colnames(sparse)
    rowfile <- paste0(stem, ".rownames")
    colfile <- paste0(stem, ".colnames")
    writeMM(sparse, stem)
    gzip(stem)
    write_lines(rows, rowfile)
    write_lines(cols, colfile)
    return(gzfile)
}
