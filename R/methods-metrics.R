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
#' @param filterCells Show only the cells that have passed filtering cutoffs.
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
        interestingGroups,
        filterCells = FALSE) {
        if (isTRUE(filterCells)) {
            object <- .applyFilterCutoffs(object)
        }
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        sampleID <- cell2sample(object)
        sampleMetadata <- sampleMetadata(object)
        colData(object) %>%
            as.data.frame() %>%
            cbind(sampleID) %>%
            left_join(sampleMetadata, by = "sampleID") %>%
            uniteInterestingGroups(interestingGroups)
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
    metadata <- sampleMetadata(object)
    suppressMessages(left_join(metrics, metadata)) %>%
        uniteInterestingGroups(interestingGroups) %>%
        column_to_rownames("cellID")
})
