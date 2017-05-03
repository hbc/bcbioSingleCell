#' Read data from bcbio-nextgen run
#'
#' @rdname read_bcbio
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run \code{bcbio-nextgen}



#' @rdname read_bcbio
#' @description Import transcript-level count data from a bcbio run into a
#'   sparse matrix
#'
#' @param strip_version Strip transcript version from identifier
#'
#' @return Sparse counts matrix
#' @export
read_bcbio_sparsecounts <- function(run, strip_version = TRUE) {
    input <- file.path(run$project_dir, "tagcounts.mtx") %>% read_sparsecounts
    if (isTRUE(strip_version)) {
        input <- strip_transcript_versions(input)
    }
    # Convert transcript-level counts to gene-level
    ensembl <- run$ensembl
    # Avoid setting rownames on data frames (tidy principle)
    indexes <- match(rownames(input), ensembl$ensembl_transcript_id)
    # collect(select(ensembl, external_gene_name))[[1]]
    gene_names <- ensembl %>%
        as.data.frame %>%
        .[indexes, "external_gene_name"]
    sparsecounts <- combine_features(input, gene_names)
    return(sparsecounts)
}
