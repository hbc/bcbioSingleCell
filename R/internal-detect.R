#' Detect scRNA-seq analysis pipeline
#'
#' Currently supports bcbio-nextgen and 10X Genomics Cell Ranger.
#'
#' @rdname detect_pipeline
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Final upload directory.
#'
#' @return Pipeline string. Stops on detection failure.
.detect_pipeline <- function(upload_dir) {
    # Search for MatrixMarket files
    matrices <- list.files(
        upload_dir, pattern = "*\\.mtx$",
        full.names = TRUE, recursive = TRUE)
    if (length(matrices) == 0) {
        stop("Failed to detect any MatrixMarket (.mtx) files")
    }
    if (any(str_detect(
        # bcbio matching by `tagcounts.mtx` in `project_dir`
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
