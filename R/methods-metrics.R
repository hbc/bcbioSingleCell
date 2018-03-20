#' Sample Barcode Metrics
#'
#' @rdname metrics
#' @name metrics
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom bcbioBase metrics
#'
#' @inheritParams general
#'
#' @param interestingGroups Interesting group, to use for colors.
#'
#' @seealso [sampleData()].
#'
#' @return `data.frame` with cellular barcodes as rows.
#'
#' @examples
#'
#' # bcbioSingleCell
#' metrics(bcb_small) %>% glimpse()
#'
#' # seurat ====
#' metrics(seurat_small) %>% glimpse()
#' metrics(pbmc_small) %>% glimpse()
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
    sampleData <- sampleData(object, interestingGroups)
    cbind(
        as.data.frame(sampleID),
        as.data.frame(colData)
    ) %>%
        rownames_to_column("cellID") %>%
        left_join(sampleData, by = "sampleID") %>%
        .sanitizeMetrics() %>%
        column_to_rownames("cellID")
}



# Methods ======================================================================
#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("bcbioSingleCell"),
    .metrics
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    function(object, interestingGroups) {
        data <- .metrics(object, interestingGroups)
        # Ensure ident column is added
        assert_are_disjoint_sets(colnames(data), "ident")
        ident <- slot(object, "ident")
        cbind(data, ident)
    }
)
