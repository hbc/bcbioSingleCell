# TODO Ensure `colData` slot gets updated upon assignment. Can we do this?



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
        # Don't allow the user to manually set sampleID column
        value[["sampleID"]] <- rownames(value)

        # Ensure the interesting groups column is not stashed
        value[["interestingGroups"]] <- NULL

        # Ensure the cell-level column data is also updated
        colData <- colData(object)
        # Require that sampleID column, defined by `cell2sample()`, is present
        assert_is_subset("sampleID", colnames(colData))
        colData <- colData[
            ,
            c("sampleID", setdiff(colnames(colData), colnames(value))),
            drop = FALSE
        ]
        colData[["cellID"]] <- rownames(colData)
        colData <- merge(
            x = colData,
            y = value,
            by = "sampleID",
            all.x = TRUE
        )
        rownames(colData) <- colData[["cellID"]]
        colData[["cellID"]] <- NULL

        # Re-slot the cell-level data
        colData <- colData[colnames(object), , drop = FALSE]
        colData(object) <- colData

        # Re-slot the sample-level data
        value[["sampleID"]] <- NULL
        metadata(object)[["sampleData"]] <- value

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
