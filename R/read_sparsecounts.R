#' Import scRNA-Seq count data
#'
#' Import transcript-level count data from a bcbio run into a sparse matrix
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run
#'
#' @return Sparse counts matrix
#' @export
read_sparsecounts <- function(run) {
    # Sparse matrix
    matfile <- file.path(run$project_dir, "tagcounts.mtx")
    if (!file.exists(matfile)) {
        stop("count matrix could not be found")
    }
    sparse <- readMM(matfile)
    # class(sparse)[1] == "dgTMatrix"

    # Rownames
    rowfile <- file.path(run$project_dir, "tagcounts.mtx.rownames")
    if (!file.exists(rowfile)) {
        stop("rownames file could not be found")
    }
    rownames(sparse) <- read_lines(rowfile)
    # Strip out transcript version numbers
    rownames(sparse) <- str_replace(rownames(sparse), "\\.\\d+", "")
    # Convert transcript-level counts to gene-level
    annotations <- ensembl_annotations(run)
    rownames(sparse) <- annotations[rownames(sparse), "external_gene_name"]

    # Colnames
    colfile <- file.path(run$project_dir, "tagcounts.mtx.colnames")
    if (!file.exists(colfile)) {
        stop("colnames file could not be found")
    }
    colnames(sparse) <- read_lines(colfile)

    # Remove rows that don't match a gene name
    sparse <- sparse[!is.na(rownames(sparse)), ]

    # Aggregate the counts by rowname
    sparse <- aggregate.Matrix(
        sparse,
        row.names(sparse),
        fun = "sum")
    # class(sparse)[1] == "dgCMatrix"

    return(as(sparse, "dgCMatrix"))
}
