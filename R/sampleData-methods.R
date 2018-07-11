#' Sample Data
#'
#' Metadata in columns describing the samples, which are defined in the
#' rownames. Similar to [colData()], which for `bcbioSingleCell` and
#' `SingleCellExperiment` objects describes cells in the columns, rather than
#' the samples.
#'
#' @name sampleData
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom basejump sampleData sampleData<-
#'
#' @inheritParams general
#' @param clean Only return factor columns.
#' @return Return either a DataFrame, data.frame, or kable.
#'
#' @examples
#' # SingleCellExperiment ====
#' x <- cellranger_small
#' sampleData(x, clean = FALSE) %>% glimpse()
#' sampleData(x, clean = TRUE) %>% glimpse()
#'
#' # Assignment support
#' sampleData(x)[["batch"]] <- 1L
#' sampleData(x) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        clean = TRUE,
        return = c("DataFrame", "data.frame", "kable")
    ) {
        object <- as(object, "SingleCellExperiment")
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        assert_is_a_bool(clean)
        return <- match.arg(return)

        data <- metadata(object)[["sampleData"]]
        if (is.null(data)) {
            return(NULL)
        }
        assert_is_data.frame(data)

        # Only return factor columns, if desired
        if (isTRUE(clean)) {
            data <- data[, vapply(data, is.factor, logical(1L)), drop = FALSE]
            # Drop remaining blacklisted columns
            setdiff <- setdiff(colnames(data), bcbioBase::metadataBlacklist)
            data <- data[, setdiff, drop = FALSE]
        } else {
            # Include `interestingGroups` column, if not NULL
            if (length(interestingGroups)) {
                data <- uniteInterestingGroups(data, interestingGroups)
            }
        }

        # Arrange rows by `sampleName` column, if defined
        if ("sampleName" %in% colnames(data)) {
            data <- data[order(data[["sampleName"]]), , drop = FALSE]
        }

        # Return
        if (return == "kable") {
            kable(as.data.frame(data), row.names = FALSE)
        } else {
            as(data, return)
        }
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("seurat"),
    getMethod("sampleData", "SingleCellExperiment")
)



# Assignment methods ===========================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "SingleCellExperiment",
        value = "DataFrame"
    ),
    function(object, value) {
        value <- sanitizeSampleData(value)
        metadata(object)[["sampleData"]] <- as.data.frame(value)
        object
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "seurat",
        value = "DataFrame"
    ),
    getMethod(
        "sampleData<-",
        signature(
            object = "SingleCellExperiment",
            value = "DataFrame"
        )
    )
)
