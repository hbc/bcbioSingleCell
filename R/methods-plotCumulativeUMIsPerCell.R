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
#' @inheritParams plotUMIsPerCell
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotCumulativeUMIsPerCell(bcb_small)
#'
#' # SingleCellExperiment ====
#' plotCumulativeUMIsPerCell(cellranger_small)
#'
#' # seurat ====
#' plotCumulativeUMIsPerCell(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotCumulativeUMIsPerCell
#' @export
setMethod(
    "plotCumulativeUMIsPerCell",
    signature("SingleCellExperiment"),
    function(
        object,
        point = c("knee", "inflection"),
        trans = c("identity", "log10", "log2", "sqrt")
    ) {
        validObject(object)
        geom <- "cumfreq"
        trans <- match.arg(trans)

        data <- metrics(object) %>%
            .[, "nUMI", drop = FALSE] %>%
            arrange(!!sym("nUMI")) %>%
            mutate(
                cumsum = cumsum(!!sym("nUMI")),
                freq = !!sym("cumsum") / max(!!sym("cumsum"))
            )

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "nUMI",
                y = "freq"
            )
        ) +
            geom_line(size = 1L) +
            scale_x_continuous(trans = trans) +
            labs(
                title = "cumulative UMIs per cell",
                x = "nUMI per cell",
                y = "cumulative frequency"
            ) +
            expand_limits(y = 0L)

        p <- .labelBarcodeRanks(
            p = p,
            object = object,
            geom = geom,
            point = point
        )

        p
    }
)



# Aliases ======================================================================
#' @rdname plotCumulativeUMIsPerCell
#' @export
plotKnee <- function(...) {
    plotCumulativeUMIsPerCell(...)  # nocov
}
