#' Read \href{https://github.com/chapmanb/bcbio-nextgen}{bcbio-nextgen} run
#' output.
#'
#' @rdname read_bcbio
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run bcbio-nextgen.



#' @rdname read_bcbio
#' @description Import transcript-level count data from a bcbio run into a
#'   sparse matrix.
#'
#' @param strip_version Strip transcript version from identifier.
#'
#' @return Sparse counts matrix.
#' @export
read_bcbio_counts <- function(run, strip_version = TRUE) {
    counts <- file.path(run$project_dir, "tagcounts.mtx") %>% read_counts
    if (isTRUE(strip_version)) {
        counts <- strip_transcript_versions(counts)
    }
    # Convert transcript-level counts to gene-level
    ensembl <- run$ensembl
    # Avoid setting rownames on data frames (tidy principle)
    indexes <- match(rownames(counts), ensembl$ensembl_transcript_id)
    # collect(select(ensembl, external_gene_name))[[1]]
    gene_names <- ensembl %>%
        as.data.frame %>%
        .[indexes, "external_gene_name"]
    counts <- aggregate_sparse_features(counts, gene_names)
    return(counts)
}
