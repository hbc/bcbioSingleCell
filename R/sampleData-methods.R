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
#' @importFrom basejump sampleData
#' @export
#'
#' @inheritParams general
#'
#' @return `DataFrame`.
#'
#' @examples
#' # SingleCellExperiment ====
#' x <- cellranger_small
#' sampleData(x) %>% glimpse()
#'
#' # Assignment support
#' sampleData(x)[["batch"]] <- 1L
#' sampleData(x) %>% glimpse()
NULL



#' @rdname sampleData
#' @name sampleData<-
#' @importFrom basejump sampleData<-
#' @export
NULL



# Methods ======================================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups
    ) {
        data <- metadata(object)[["sampleData"]]
        if (is.null(data)) {
            return(NULL)
        }
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        if (length(interestingGroups)) {
            data <- uniteInterestingGroups(data, interestingGroups)
        }
        as(data, "DataFrame")
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
        value[["interestingGroups"]] <- NULL
        metadata(object)[["sampleData"]] <- value
        object
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "SingleCellExperiment",
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
