#' Detect scRNA-Seq Analysis Pipeline
#'
#' Currently supports bcbio-nextgen and 10X Genomics Cell Ranger.
#'
#' @rdname detectPipeline-internal
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param uploadDir Final upload directory.
#'
#' @return Pipeline string. Stops on detection failure.
.detectPipeline <- function(uploadDir) {
    # Search for MatrixMarket files
    matrices <- list.files(
        uploadDir, pattern = "*\\.mtx$",
        full.names = TRUE, recursive = TRUE)
    if (length(matrices) == 0L) {
        stop("Failed to detect any MatrixMarket (.mtx) files")
    }
    if (any(str_detect(
        # bcbio matching by `tagcounts.mtx` in `projectDir`
        matrices, "\\d{4}-\\d{2}-\\d{2}_[^/]+/tagcounts\\.mtx$"))) {
        "bcbio"
    } else if (any(str_detect(
        # cellranger matching by raw `matrix.mtx`
        matrices, "outs/raw_gene_bc_matrices/[^/]+/matrix\\.mtx$"))) {
        "cellranger"
    } else {
        stop("Pipeline detection failed")
    }
}
