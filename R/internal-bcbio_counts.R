#' Import transcript-level count data from a bcbio run into a sparse matrix
#'
#' @rdname counts
#' @author Michael Steinbaugh
#'
#' @param strip_version Strip transcript version from identifier.
#'
#' @return Sparse counts matrix.
.bcbio_counts <- function(project_dir, strip_version = TRUE) {
    counts <- file.path(project_dir, "tagcounts.mtx") %>% .sparse_counts
    if (isTRUE(strip_version)) {
        counts <- .strip_transcript_versions(counts)
    }
    # Convert transcript-level counts to gene-level
    # FIXME Need to switch to tx2gene method (annotables?)
    stop("Ensembl method not fixed yet")
    indexes <- match(rownames(counts), ensembl[["ensembl_transcript_id"]])
    gene_names <- ensembl %>%
        as.data.frame %>%
        .[indexes, "external_gene_name"]
    .aggregate_sparse_features(counts, gene_names)
}
