#' Metrics per Sample
#'
#' Calculate summary statistics per sample.
#'
#' @name metricsPerSample
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param f Mathematical function to apply. Defaults to "`sum`".
#'
#' @return `data.frame`.
#'
#' @examples
#' # SingleCellExperiment ====
#' # mean
#' x <- metricsPerSample(bcb_small, f = "mean")
#' head(x)
#'
#' # median
#' x <- metricsPerSample(bcb_small, f = "median")
#' head(x)
#'
#' # sum
#' x <- metricsPerSample(bcb_small, f = "sum")
#' head(x)
NULL



# Methods ======================================================================
#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("SingleCellExperiment"),
    function(
        object,
        f = c("mean", "median", "sum")
    ) {
        f <- match.arg(f)
        inform(paste("Calculating", f, "per sample"))
        fxn <- get(f)
        assert_is_function(fxn)
        assert_is_subset("sampleName", colnames(metrics(object)))
        metrics <- metrics(object)
        if (f == "sum") {
            # Sum only the `n*` columns containing counts
            data <- select(metrics, matches("^n[A-Z]"))
        } else {
            # Summarize all numeric columns
            data <- select_if(metrics, is.numeric)
        }
        assert_is_non_empty(data)
        data[["rowname"]] <- metrics[["sampleName"]]
        data %>%
            group_by(!!sym("rowname")) %>%
            summarize_all(fxn) %>%
            arrange(!!sym("rowname")) %>%
            ungroup() %>%
            as.data.frame() %>%
            column_to_rownames()
    }
)



#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("seurat"),
    getMethod("metricsPerSample", "SingleCellExperiment")
)
