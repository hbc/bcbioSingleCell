#' Read data from bcbio-nextgen run
#'
#' @rdname read_bcbio
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run \code{bcbio-nextgen}
#' @param save Save data



#' @rdname read_bcbio
#' @description Import metadata
#'
#' @return Metadata data frame
#' @export
read_bcbio_metadata <- function(
    run,
    save = FALSE) {
    check_run(run)
    metadata <- list.files(run$config_dir,
                           pattern = ".csv$",
                           full.names = TRUE) %>%
        read_csv(col_types = cols()) %>%
        set_names_snake %>%
        arrange_(.dots = "description") %>%
        as.data.frame %>%
        set_rownames(.$description)

    if (isTRUE(save)) {
        save(metadata, file = "data/metadata.rda")
        write_csv(metadata, "meta/metadata.csv")
    }

    return(metadata)
}



#' Import single-cell RNA-seq count data
#'
#' Import transcript-level count data from a bcbio run into a sparse matrix
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run \code{bcbio-nextgen} run
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
    annotations <- ensembl_annotations(run)
    # Avoid setting rownames on data frames (tidy principle)
    indexes <- match(rownames(input), annotations$ensembl_transcript_id)
    # collect(select(annotations, external_gene_name))[[1]]
    gene_names <- annotations %>%
        as.data.frame %>%
        .[indexes, "external_gene_name"]
    sparsecounts <- combine_features(input, gene_names)
    return(sparsecounts)
}
