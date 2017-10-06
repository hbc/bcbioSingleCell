#' Plot Percentage of Zeros vs. Library Depth
#'
#' @rdname plotZerosVsDepth
#' @name plotZerosVsDepth
#' @family Quality Control Metrics
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [ggplot].
NULL



# Constructors ====
.plotZerosVsDepth <- function(object) {
    counts <- assay(object)
    metrics <- metrics(object)
    # using a logical matrix is much faster than a == 0 comparison
    present <- as(counts, "lgCMatrix")

    # The `counts == 0` operation here blows up memory on large datasets
    df <- data.frame(
        dropout = (nrow(present) - Matrix::colSums(present)) / nrow(present),
        depth = Matrix::colSums(counts),
        description = metrics[["description"]])
    p <- ggplot(
        df,
        mapping = aes_(
            x = ~depth,
            y = ~dropout * 100)
    ) +
        geom_point(size = 0.8, alpha = 0.3) +
        scale_x_log10() +
        labs(x = "library size",
             y = "% genes zero")
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p +
            facet_wrap(facets = "description")
    }
    p
}



# Methods ====
#' @rdname plotZerosVsDepth
#' @export
setMethod("plotZerosVsDepth", "bcbioSingleCellANY", .plotZerosVsDepth)
