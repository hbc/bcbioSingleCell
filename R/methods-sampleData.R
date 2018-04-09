#' Sample Data
#'
#' Metadata in columns describing the samples, which are defined in the
#' rownames. Similar to [colData()], which for `bcbioSingleCell` and
#' `SingleCellExperiment` objects describes cells in the columns, rather than
#' the samples.
#'
#' @name sampleData
#' @family S4 Class Definitions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase sampleData sampleData<-
#'
#' @inheritParams general
#'
#' @return `data.frame`.
#'
#' @examples
#' # bcbioSingleCell ====
#' x <- bcb_small
#' sampleData(x) %>% glimpse()
#' sampleData(x)[["batch"]] <- 1L
#' sampleData(x) %>% glimpse()
#'
#' # seurat ====
#' x <- seurat_small
#' sampleData(x) %>% glimpse()
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
        return = c("DataFrame", "data.frame", "kable")
    ) {
        validObject(object)
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        return <- match.arg(return)
        data <- metadata(object)[["sampleData"]]
        assert_is_any_of(data, c("DataFrame", "data.frame"))
        data <- uniteInterestingGroups(data, interestingGroups)
        data <- sanitizeSampleData(data)
        assertHasRownames(data)
        if (return == "kable") {
            blacklist <- c("description", "fileName", "sampleID")
            data %>%
                as.data.frame() %>%
                .[, setdiff(colnames(.), blacklist), drop = FALSE] %>%
                # Ensure `sampleName` is first
                .[, unique(c("sampleName", colnames(.))), drop = FALSE] %>%
                # Arrange by `sampleName`
                .[order(.[["sampleName"]]), , drop = FALSE] %>%
                kable(row.names = FALSE)
        } else {
            as(data, return)
        }
    }
)



# Assignment methods ===========================================================
#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "SingleCellExperiment",
        value = "ANY"
    ),
    function(object, value) {
        # TODO Switch to storing as DataFrame
        value <- as.data.frame(value)
        # Ensure all columns are factors
        value <- sanitizeSampleData(value)
        metadata(object)[["sampleData"]] <- value
        object
    }
)
