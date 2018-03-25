#' Column Data
#'
#' @name colData
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom SummarizedExperiment colData
#'
#' @inheritParams general
#' @param return Return data as `data.frame` or `DataFrame`.
#'
#' @return Data describing the columns (cells).
#'
#' @examples
#' # seurat ====
#' colData(pbmc_small)
NULL



# Methods ======================================================================
#' @rdname colData
#' @export
setMethod(
    "colData",
    signature("bcbioSingleCell"),
    function(x, return = c("data.frame", "DataFrame")) {
        return <- match.arg(return)
        data <- slot(x, "colData")
        as(data, return)
    }
)



#' @rdname colData
#' @export
setMethod(
    "colData",
    signature("seurat"),
    function(x, return = c("data.frame", "DataFrame")) {
        return <- match.arg(return)

        data <- slot(x, "meta.data")
        assert_is_data.frame(data)
        # Sanitize column names
        data <- camel(data)

        # Legacy: ensure colData doesn't contain sample metadata
        meta <- sampleData(x)
        drop <- intersect(
            x = colnames(data),
            y = c(colnames(meta), "cellularBarcode")
        )
        if (length(drop)) {
            data <- data[, setdiff(colnames(data), drop)]
        }

        data %>%
            rownames_to_column() %>%
            mutate_if(is.character, as.factor) %>%
            mutate_if(is.factor, droplevels) %>%
            column_to_rownames() %>%
            as(return)
    }
)
