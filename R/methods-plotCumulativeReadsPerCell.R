#' Plot Cumulative Frequency of Reads Per Cell
#'
#' Cumulative frequency plot that helps visualize the inflection point that is
#' present in single-cell barcode distributions.
#'
#' @name plotCumulativeReadsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams inflectionPoint
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotCumulativeReadsPerCell(bcb_small)
#'
#' # seurat ====
#' plotCumulativeReadsPerCell(pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotCumulativeReadsPerCell
#' @export
setMethod(
    "plotCumulativeReadsPerCell",
    signature("bcbioSingleCell"),
    function(
        object,
        maxCells = 100000L
    ) {
        assertIsAnImplicitInteger(maxCells)

        metrics <- metrics(object)
        countsCol <- c("nCount", "nUMI")
        assert_are_intersecting_sets(countsCol, colnames(metrics))
        countsCol <- intersect(countsCol, colnames(metrics))[[1L]]

        data <- metrics %>%
            as_tibble() %>%
            rownames_to_column() %>%
            .[, c("rowname", countsCol)] %>%
            arrange(!!sym(countsCol)) %>%
            mutate(
                cumsum = cumsum(!!sym(countsCol)),
                freq = !!sym("cumsum") / max(!!sym("cumsum"))
            )
        inflection <- inflectionPoint(object, maxCells = maxCells)

        ggplot(
            data = data,
            mapping = aes_string(
                x = countsCol,
                y = "freq"
            )
        ) +
            geom_line() +
            .qcCutoffLine(xintercept = inflection) +
            labs(
                y = "cumulative frequency"
            ) +
            expand_limits(x = 0L, y = 0L)
    }
)



# Aliases ======================================================================
#' @rdname plotCumulativeReadsPerCell
#' @export
plotKnee <- function(...) {
    plotCumulativeReadsPerCell(...)  # nocov
}
