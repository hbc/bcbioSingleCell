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
#' plotCumulativeReadsPerCell(cellranger_small)
#'
#' # seurat ====
#' plotCumulativeReadsPerCell(seurat_small)
#' plotCumulativeReadsPerCell(pbmc_small)
NULL



# Constructors =================================================================
.plotCumulativeReadsPerCell <- function(
    object,
    maxCells = 100000L
) {
    validObject(object)
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
            title = "cumulative reads per cell",
            subtitle = paste("inflection", inflection, sep = " = "),
            y = "cumulative frequency"
        ) +
        expand_limits(y = 0L)
}


# Methods ======================================================================
#' @rdname plotCumulativeReadsPerCell
#' @export
setMethod(
    "plotCumulativeReadsPerCell",
    signature("bcbioSingleCell"),
    .plotCumulativeReadsPerCell
)



#' @rdname plotCumulativeReadsPerCell
#' @export
setMethod(
    "plotCumulativeReadsPerCell",
    signature("seurat"),
    .plotCumulativeReadsPerCell
)



# Aliases ======================================================================
#' @rdname plotCumulativeReadsPerCell
#' @export
plotKnee <- function(...) {
    plotCumulativeReadsPerCell(...)  # nocov
}
