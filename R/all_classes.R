#' Class that contains bcbio run information
#'
#' `bcbioSCDataSet` is a subclass of [SummarizedExperiment] designed to store a
#' single-cell RNA-seq analysis. This class contains experiment metadata, raw
#' counts, normalilzed counts, and summary statistics for each sample analyzed.
#'
#' Methods for this objects ...
#'
#' `metadata` contains ...
#'
#' @rdname bcbioSCDataSet
#' @keywords internal
#' @aliases bcbioSCDataSet-class
#' @export
bcbioSCDataSet <- setClass(
    "bcbioSCDataSet",
    contains = "SummarizedExperiment",
    slots = c("sparse" = "dgCMatrix"))
# setValidity("bcbioSCDataSet", function(object) TRUE)

# Store the metrics in the assays slot!
# https://github.com/satijalab/seurat/blob/master/R/seurat.R

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
#' We recommend loading the bcbio-nextgen run as a remote connection over
#' `sshfs`.
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
