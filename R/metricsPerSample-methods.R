#' Metrics per Sample
#'
#' Calculate summary statistics per sample.
#'
#' @name metricsPerSample
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param f `string`. Mathematical function name to apply. Defaults to "`sum`"
#'   but "`mean`" and "`median`" are also supported.
#'
#' @return `data.frame`.
#'
#' @examples
#' # mean
#' x <- metricsPerSample(indrops_small, f = "mean")
#' head(x)
#'
#' # median
#' x <- metricsPerSample(indrops_small, f = "median")
#' head(x)
#'
#' # sum
#' x <- metricsPerSample(indrops_small, f = "sum")
#' head(x)
NULL



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
        message(paste("Calculating", f, "per sample"))
        fxn <- get(f)
        assert_is_function(fxn)
        metrics <- metrics(object)
        assert_is_subset("sampleName", colnames(metrics))
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
