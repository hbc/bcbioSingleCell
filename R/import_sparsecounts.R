#' Import scRNA-Seq count data
#'
#' Import transcript-level count data from a bcbio run into a sparse matrix
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import readr
#' @importFrom Matrix readMM
#' @importFrom Matrix.utils aggregate.Matrix
#'
#' @param bcbio bcbio run object
#' @param annotations Ensembl annotations data frame
#' @param tx2gene Convert transcript-level counts to gene-level
#' @param strip_versions Strip transcript versions from the matrix
#'   before aggregating. Only applies to transcript-level output.
#' @param print Whether to print \code{Dimnames}
#'
#' @return Sparse counts matrix
#' @export
#'
#' @examples
#' \dontrun{
#' import_sparsecounts(bcbio)
#' }
import_sparsecounts <- function(
    bcbio,
    annotations,
    tx2gene = TRUE,
    strip_versions = FALSE,
    print = FALSE) {
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

    # Strip transcript versions for better aggregation and/or gene matching
    if (isTRUE(strip_versions)) {
        sparse <- strip_transcript_versions(sparse)
    }

    # Convert transcript-level counts to gene-level
    if (isTRUE(tx2gene)) {
        tx2gene <- annotations %>%
            dplyr::select_(.dots = c("ensembl_transcript_id",
                                     "ensembl_gene_id")) %>%
            dplyr::distinct(.) %>%
            set_rownames("ensembl_transcript_id")
        rownames(sparse) <- tx2gene[rownames(sparse), "ensembl_gene_id"]
    }

    # Remove rows that don't match to a gene identifier
    sparse <- sparse[!is.na(rownames(sparse)), ]

    # Aggregate the counts by rowname
    sparse <- Matrix.utils::aggregate.Matrix(
        sparse,
        row.names(sparse),
        fun = "sum"
    )

    # Diagnostic output of row and column names
    if (isTRUE(print)) {
        print(sparse@Dimnames[[1]][1:3])
        print(sparse@Dimnames[[2]][1:3])
    }

    return(sparse)
}
