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

##' Convert a transcript matrix to a gene symbol matrix
##'
##' @param matrix a matrix of counts
##' @param fromto a dataframe where "from" is the original id and "to" is the id to convert to
##' @param strip strip transcript versions from the matrix before aggregating
##' @return a matrix of counts aggregated by summing the "to" ids
##' @import dplyr
##' @importFrom Matrix.utils aggregate.Matrix
##' @author Rory Kirchner
##' @export
tx2symbol <- function(matrix, fromto, strip = FALSE) {
    if (strip) {
        matrix <- strip_transcript_versions(matrix)
    }
    lookup <- data.frame(from = rownames(matrix)) %>%
        dplyr::left_join(fromto, by = "from")
    rownames(matrix) <- lookup$to
    matrix <- matrix[!is.na(rownames(matrix)), ]
    matrix <- Matrix.utils::aggregate.Matrix(matrix,
                                             row.names(matrix),
                                             fun = "sum")
    return(matrix)
}

##' Strip transcript versions off of the rownames of a dataframe or matrix
##'
##' @param matrix a matrix or dataframe
##' @return a matrix or dataframe with the transcript versions stripped off
##' @importFrom tools file_path_sans_ext
##' @author Rory Kirchner
##' @export
strip_transcript_versions <- function(matrix) {
    transcripts <- rownames(matrix)
    transcripts <- tools::file_path_sans_ext(transcripts)
    rownames(matrix) <- transcripts
    return(matrix)
}
