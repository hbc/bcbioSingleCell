#' Column Data
#'
#' @family S4 Class Definition
#' @author Michael Steinbaugh
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
        sampleID <- DataFrame("sampleID" = cell2sample(x))
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
        data
    }
)
