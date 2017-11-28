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
#' @importFrom Matrix colSums
.plotZerosVsDepth <- function(object) {
    counts <- assay(object)
    metrics <- metrics(object)
    # Using a logical matrix is much faster and more memory efficient than
    # `counts == 0` comparison
    present <- as(counts, "lgCMatrix")
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

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (!is.null(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}



# Methods ====
#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("bcbioSingleCell"),
    .plotZerosVsDepth)
