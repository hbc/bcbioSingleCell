#' Sample Barcode Metrics
#'
#' @rdname metrics
#' @name metrics
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump metrics
#'
#' @inheritParams AllGenerics
#'
#' @param interestingGroups Interesting group, to use for colors.
#'
#' @seealso [sampleMetadata()].
#'
#' @return [data.frame] with cellular barcodes as rows.
NULL



# Methods ====
#' @rdname metrics
#' @importFrom basejump uniteInterestingGroups
#' @importFrom dplyr left_join mutate_all
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
setMethod(
    "metrics",
    signature("bcbioSingleCell"),
    function(
        object,
        interestingGroups) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }

        colData <- colData(object)
        # Check for malformed colData. Safe to remove in a future update.
        if ("sampleID" %in% colnames(colData)) {
            warning("Sample ID shouldn't be set in colData",
                 call. = FALSE)
            colData[["sampleID"]] <- NULL
        }

        sampleID <- cell2sample(object)
        sampleMetadata <- sampleMetadata(
            object,
            interestingGroups = interestingGroups)

        cbind(
            as.data.frame(sampleID),
            as.data.frame(colData)
        ) %>%
            left_join(sampleMetadata, by = "sampleID")
    })



#' @rdname metrics
#' @importFrom basejump camel
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
        interestingGroups <- basejump::interestingGroups(object)
    }
    metrics <- slot(object, "meta.data") %>%
        as.data.frame() %>%
        camel(strict = FALSE) %>%
        rownames_to_column("cellID") %>%
        mutate_if(is.character, as.factor)
    metrics <- suppressMessages(left_join(
        metrics,
        sampleMetadata(object, interestingGroups = interestingGroups)
    ))
    metrics %>%
        uniteInterestingGroups(interestingGroups) %>%
        column_to_rownames("cellID")
})
