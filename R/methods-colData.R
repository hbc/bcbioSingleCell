#' Column Data
#'
#' @rdname colData
#' @name colData
#' @author Michael Steinbaugh
#'
#' @importFrom SummarizedExperiment colData
#'
#' @inherit general
NULL



# Methods ======================================================================
#' @rdname colData
#' @importFrom dplyr mutate_if
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
setMethod(
    "colData",
    signature("seurat"),
    function(x) {
        data <- slot(x, "meta.data") %>%
            as.data.frame() %>%
            # Ensure the column names get sanitized
            camel()

        # Legacy: ensure colData doesn't contain sample metadata
        meta <- sampleData(x)
        drop <- intersect(
            x = colnames(data),
            y = c(
                colnames(meta),
                "cellularBarcode"
            )
        )
        if (length(drop)) {
            data <- data %>%
                .[, setdiff(colnames(.), drop)]
        }

        data %>%
            rownames_to_column() %>%
            mutate_if(is.character, as.factor) %>%
            mutate_if(is.factor, droplevels) %>%
            column_to_rownames()
    })
