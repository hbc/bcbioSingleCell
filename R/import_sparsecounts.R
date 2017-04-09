#' Import scRNA-Seq count data
#'
#' Import transcript-level count data from a bcbio run into a sparse matrix
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import readr
#' @import Matrix
#' @importFrom stringr str_replace
#' @importFrom Matrix.utils aggregate.Matrix
#'
#' @param bcbio bcbio run object
#' @param annotations Ensembl annotations data frame
#'
#' @return Sparse counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' import_sparsecounts(bcbio)
#' }
import_sparsecounts <- function(bcbio, annotations) {
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
    # strip out ensembl transcript version numbers
    rownames(sparse) <- stringr::str_replace(rownames(sparse), "\\.\\d+", "")

    colnames(sparse) <- readr::read_lines(colfile)

    # Convert transcript-level counts to gene-level
    rownames(sparse) <- annotations[rownames(sparse), "gene_name"]

    # Remove rows that don't match a gene name
    sparse <- sparse[!is.na(rownames(sparse)), ]

    # Aggregate the counts by rowname
    sparse <- Matrix.utils::aggregate.Matrix(
        sparse,
        row.names(sparse),
        fun = "sum"
    )
    return(as(sparse, "dgCMatrix"))
}
