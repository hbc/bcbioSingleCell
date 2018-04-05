#' Plot Cumulative Frequency of UMIs Per Cell
#'
#' Cumulative frequency plot that helps visualize the inflection point that is
#' present in single-cell barcode distributions.
#'
#' These calculates apply to all cellular barcodes in the dataset and are not
#' grouped by sample.
#'
#' @name plotCumulativeUMIsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' # SingleCellExperiment ====
#' plotCumulativeUMIsPerCell(bcb_small)
#' plotCumulativeUMIsPerCell(cellranger_small)
NULL



# Methods ======================================================================
#' @rdname plotCumulativeUMIsPerCell
#' @export
setMethod(
    "plotCumulativeUMIsPerCell",
    signature("SingleCellExperiment"),
    function(
        object,
        trans = c("identity", "log10", "log2", "sqrt")
    ) {
        validObject(object)
        trans <- match.arg(trans)

        data <- metrics(object) %>%
            .[, "nUMI", drop = FALSE] %>%
            arrange(!!sym("nUMI")) %>%
            mutate(
                cumsum = cumsum(!!sym("nUMI")),
                freq = !!sym("cumsum") / max(!!sym("cumsum"))
            )

        inflection <- inflectionPoint(object)

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "nUMI",
                y = "freq"
            )
        ) +
            geom_line() +
            scale_x_continuous(trans = trans) +
            labs(
                title = "cumulative UMIs per cell",
                subtitle = paste("inflection", inflection, sep = " = "),
                x = "nUMI per cell",
                y = "cumulative frequency"
            ) +
            expand_limits(y = 0L)

        if (inflection > 0L) {
            p <- p +
                .qcCutoffLine(
                    xintercept = inflection,
                    color = inflectionColor
                )
        }

        p
    }
)



# Aliases ======================================================================
#' @rdname plotCumulativeUMIsPerCell
#' @export
plotKnee <- function(...) {
    plotCumulativeUMIsPerCell(...)  # nocov
}
