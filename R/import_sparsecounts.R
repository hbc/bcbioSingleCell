#' Import scRNA-Seq count data
#'
#' Import transcript-level count data from a bcbio run into a sparse matrix
#'
#' @author Michael Steinbaugh
#'
#' @keywords internal
#'
#' @import readr
#' @importFrom Matrix readMM
#'
#' @param bcbio bcbio run object
#' @param print Whether to print \code{Dimnames}
#'
#' @return Sparse counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' import_sparsecounts(bcbio)
#' }
import_sparsecounts <- function(bcbio, print = TRUE) {
    matfile <- file.path(bcbio$project_dir, "tagcounts.mtx")
    if (!file.exists(matfile)) {
        stop("Count matrix could not be found.")
    }
    rowfile <- file.path(bcbio$project_dir, "tagcounts.mtx.rownames")
    if (!file.exists(rowfile)) {
        stop("Row names file could not be found.")
    }
    colfile <- file.path(bcbio$project_dir, "tagcounts.mtx.colnames")
    if (!file.exists(colfile)) {
        stop("Column names file could not be found.")
    }

    sparse <- Matrix::readMM(matfile)
    rownames(sparse) <- readr::read_lines(rowfile)
    colnames(sparse) <- readr::read_lines(colfile)

    if (isTRUE(print)) {
        print(counts@Dimnames[[1]][1:3])
        print(counts@Dimnames[[2]][1:3])
    }

    return(sparse)
}
