# Constructor
.bcbioSCDataSet <- function(se, sparse) {
    se <- as(se, "SummarizedExperiment")
    if (!is(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object")
    }
    object <- new("bcbioSCDataSet", se)
    object@sparse <- sparse
    # bcbio(object, "sparse") <- sparse
    object
}



#' Load bcbio-nextgen run
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh
#'
#' @param upload_dir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param metadata_file Sample barcode metadata file.
#' @param organism Organism name, following Ensembl/Biomart conventions. Must be
#'   lowercase and one word (e.g. hsapiens). This will be detected automatically
#'   for common reference genomes.
#' @param intgroup Character vector of interesting groups. First entry is used
#'   for plot colors during quality control (QC) analysis. Entire vector is used
#'   for PCA and heatmap QC functions.
#'
#' @return [bcbioSCDataSet].
#' @export
load_run <- function(
    upload_dir = "final",
    # [fix] External metadata file is required until YAML support added
    metadata_file,
    organism,
    intgroup = "sample_name") {
    run <- load_run_as_list(
        upload_dir,
        metadata_file = metadata_file,
        organism = organism,
        intgroup = intgroup)

    # Use the metrics for colData
    counts <- run$sparse
    colData <- run$metrics[colnames(counts), ]
    identical(colnames(run$sparse), rownames(colData))

    run$counts <- NULL
    run$sparse <- NULL

    se <- SummarizedExperiment(
        assays = SimpleList(counts = counts),
        colData = colData,
        metadata = run)

    # Print the object size (for testing)
    object.size(se) %>% print(units = "auto")

    .bcbioSCDataSet(se, run)
}
