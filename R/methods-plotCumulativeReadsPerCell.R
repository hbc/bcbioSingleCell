#' Plot Cumulative Frequency of Reads Per Cell
#'
#' Cumulative frequency plot that helps visualize the inflection point that is
#' present in single-cell barcode distributions.
#'
#' These calculates apply to all cellular barcodes in the dataset and are not
#' grouped by sample.
#'
#' @name plotCumulativeReadsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
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
    trans = c("log10", "log2", "sqrt", "identity")
) {
    validObject(object)
    trans <- match.arg(trans)

    data <- .readsPerCell(object) %>%
        ungroup() %>%
        .[, "nCount", drop = FALSE] %>%
        arrange(!!sym("nCount")) %>%
        mutate(
            cumsum = cumsum(!!sym("nCount")),
            freq = !!sym("cumsum") / max(!!sym("cumsum"))
        )

    cutoff <- metadata(object)[["cellularBarcodeCutoff"]]
    if (!is.numeric(cutoff)) {
        cutoff <- 0L
    }
    inflection <- inflectionPoint(object)

    p <- ggplot(
        data = data,
        mapping = aes_string(
            x = "nCount",
            y = "freq"
        )
    ) +
        geom_line() +
        scale_x_continuous(trans = trans) +
        labs(
            title = "cumulative reads per cell",
            subtitle = paste(
                paste("cutoff", cutoff, sep = " = "),
                paste("inflection", inflection, sep = " = "),
                sep = "\n"
            ),
            x = "reads per cell",
            y = "cumulative frequency"
        ) +
        expand_limits(y = 0L)

    if (inflection > 0L) {
        p <- p + .qcCutoffLine(xintercept = inflection, color = inflectionColor)
    }

    p
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
