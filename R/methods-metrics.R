#' Sample Barcode Metrics
#'
#' @name metrics
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom bcbioBase metrics
#'
#' @inheritParams general
#'
#' @seealso [sampleData()].
#'
#' @return `data.frame` with cellular barcodes as rows.
#'
#' @examples
#' # bcbioSingleCell ====
#' metrics(bcb_small) %>% glimpse()
#'
#' # seurat ====
#' metrics(seurat_small) %>% glimpse()
#' metrics(pbmc_small) %>% glimpse()
NULL



# Constructors =================================================================
.metrics <- function(object, interestingGroups) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }

    colData <- colData(object, return = "data.frame")
    # Bind `sampleID` column to colData
    sampleID <- cell2sample(object)
    assert_are_identical(rownames(colData), names(sampleID))
    colData <- cbind(colData, sampleID)
    colData <- rownames_to_column(colData, "cellID")

    sampleData <- sampleData(object)

    merge(
        x = colData,
        y = sampleData,
        by = "sampleID",
        all.x = TRUE
    ) %>%
        .sanitizeMetrics() %>%
        # Ensure the metrics (colData) columns appear first
        .[, unique(c(colnames(colData), colnames(.)))] %>%
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
