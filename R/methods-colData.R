#' Column Data
#'
#' @name colData
#' @family S4 Class Definition
#' @author Michael Steinbaugh
#'
#' @importFrom SummarizedExperiment colData
#'
#' @inheritParams general
#'
#' @examples
#' # bcbioSingleCell ====
#' colData(bcb_small) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname colData
#' @export
setMethod(
    "colData",
    signature("bcbioSingleCell"),
    function(x) {
        colData <- slot(x, "colData")
        cell2sample <- metadata(x)[["cell2sample"]]
        assert_is_factor(cell2sample)
        sampleID <- DataFrame("sampleID" = cell2sample)
        colData <- cbind(colData, sampleID)
        colData[["rowname"]] <- rownames(colData)
        sampleData <- as(metadata(x)[["sampleData"]], "DataFrame")
        data <- merge(
            x = colData,
            y = sampleData,
            by = "sampleID",
            all.x = TRUE
        )
        # Ensure the numeric metrics columns appear first
        data <- data[, unique(c(colnames(colData), colnames(data)))]
        rownames(data) <- data[["rowname"]]
        data[["rowname"]] <- NULL
        # Add `interestingGroups` column
        interestingGroups <- metadata(x)[["interestingGroups"]]
        data <- uniteInterestingGroups(data, interestingGroups)
        data
    }
)
