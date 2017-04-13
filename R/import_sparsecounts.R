#' Import scRNA-Seq count data
#'
#' Import transcript-level count data from a bcbio run into a sparse matrix
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
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
    sparse <- readMM(matfile)
    rownames(sparse) <- read_lines(rowfile)
    # strip out ensembl transcript version numbers
    rownames(sparse) <- str_replace(rownames(sparse), "\\.\\d+", "")

    colnames(sparse) <- read_lines(colfile)

    # Convert transcript-level counts to gene-level
    rownames(sparse) <- annotations[rownames(sparse), "gene_name"]

    # Remove rows that don't match a gene name
    sparse <- sparse[!is.na(rownames(sparse)), ]

    # Aggregate the counts by rowname
    sparse <- aggregate.Matrix(
        sparse,
        row.names(sparse),
        fun = "sum"
    )
    return(as(sparse, "dgCMatrix"))
}
