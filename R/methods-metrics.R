#' Sample Barcode Metrics
#'
#' @rdname metrics
#' @name metrics
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom bcbioBase metrics
#'
#' @inheritParams AllGenerics
#'
#' @param interestingGroups Interesting group, to use for colors.
#'
#' @seealso [sampleMetadata()].
#'
#' @return [data.frame] with cellular barcodes as rows.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' metrics(bcb) %>% glimpse()
#'
#' # seurat
#' metrics(seurat) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname metrics
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom dplyr left_join
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
setMethod(
    "metrics",
    signature("bcbioSingleCell"),
    function(
        object,
        interestingGroups) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        colData <- colData(object)
        sampleID <- cell2sample(object)
        # Check for cell ID dimension mismatch. This can happen if `cell2sample`
        # mapping isn't updated inside the metadata slot.
        if (!identical(rownames(colData), names(sampleID))) {
            stop("'cellID' mismatch between 'colData' and 'cell2sample'",
                 call. = FALSE)
        }
        sampleMetadata <- sampleMetadata(
            object,
            interestingGroups = interestingGroups)
        metrics <- cbind(
            as.data.frame(sampleID),
            as.data.frame(colData)
        ) %>%
            rownames_to_column("cellID") %>%
            left_join(sampleMetadata, by = "sampleID") %>%
            column_to_rownames("cellID")
        metrics
    })



#' @rdname metrics
#' @importFrom bcbioBase camel
#' @importFrom dplyr mutate_if
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    function(
        object,
        interestingGroups) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    metrics <- slot(object, "meta.data") %>%
        as.data.frame() %>%
        camel(strict = FALSE) %>%
        rownames_to_column("cellID") %>%
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels)
    # Join columns can vary here, so suppress message
    metrics <- suppressMessages(left_join(
        metrics,
        sampleMetadata(object, interestingGroups = interestingGroups)
    )) %>%
        column_to_rownames("cellID") %>%
        uniteInterestingGroups(interestingGroups)
    metrics
})
