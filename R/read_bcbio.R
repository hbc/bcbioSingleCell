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
read_metadata <- function(
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
read_bcbio_sparsecounts <- function(run) {
    # Sparse matrix
    matfile <- file.path(run$project_dir, "tagcounts.mtx")
    sparse <- read_sparsecounts(matfile)

    # Strip out transcript version numbers
    rownames(sparse) <- str_replace(rownames(sparse), "\\.\\d+", "")

    # Convert transcript-level counts to gene-level
    annotations <- ensembl_annotations(run)
    geneids <- annotations[rownames(sparse), "external_gene_name"]
    sparse <- combine_features(sparse, geneids)
    return(sparse)
}
