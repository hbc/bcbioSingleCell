#' Sample Barcode Metrics
#'
#' @rdname metrics
#' @name metrics
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom bcbioBase metrics
#'
#' @inheritParams general
#'
#' @param interestingGroups Interesting group, to use for colors.
#'
#' @seealso [sampleMetadata()].
#'
#' @return [data.frame] with cellular barcodes as rows.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' metrics(bcb) %>% glimpse()
#'
#' # seurat
#' metrics(seurat) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom bcbioBase uniteInterestingGroups
#' @importFrom dplyr left_join
#' @importFrom tibble column_to_rownames rownames_to_column
.metrics <- function(object, interestingGroups) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    colData <- colData(object)
    sampleID <- cell2sample(object)
    assert_are_identical(rownames(colData), names(sampleID))
    sampleMetadata <- sampleMetadata(
        object,
        interestingGroups = interestingGroups)
    cbind(
        as.data.frame(sampleID),
        as.data.frame(colData)
    ) %>%
        # Need to use `cellID` column for rownames, because that is the
        # column set in `sampleMetadata`
        rownames_to_column("cellID") %>%
        left_join(sampleMetadata, by = "sampleID") %>%
        .sanitizeMetrics() %>%
        column_to_rownames("cellID")
}



# Methods ======================================================================
#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("bcbioSingleCell"),
    .metrics)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    .metrics)
