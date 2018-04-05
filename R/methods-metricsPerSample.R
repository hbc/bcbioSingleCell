#' Metrics per Sample
#'
#' @name metricsPerSample
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param f Mathematical function to apply. Defaults to "`sum`".
#'
#' @return `tbl_df` containing summary statistics.
#'
#' @examples
#' # bcbioSingleCell ====
#' metricsPerSample(bcb_small)
#'
#' # seurat ====
#' metricsPerSample(seurat_small)
NULL



# Methods ======================================================================
#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("SingleCellExperiment"),
    function(
        object,
        f = c("sum", "mean", "median")
    ) {
        f <- match.arg(f)
        inform(paste("Calculating", f, "per sample"))
        f <- get(f)
        assert_is_function(f)
        assert_is_subset("sampleName", colnames(metrics(object)))
        metrics <- metrics(object)
        # Summarize only the `n*` columns containing counts
        nCols <- grep("^n[A-Z]", colnames(metrics), value = TRUE)
        assert_is_non_empty(nCols)
        data <- metrics[, c("sampleName", nCols)]
        data %>%
            group_by(!!sym("sampleName")) %>%
            summarize_all(f)
    }
)



#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("seurat"),
    getMethod("metricsPerSample", "SingleCellExperiment")
)
