#' Plot Percentage of Zeros vs. Library Depth
#'
#' @rdname plotZerosVsDepth
#' @name plotZerosVsDepth
#' @family Quality Control Metrics
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @return [ggplot].
NULL



# Constructors ====
.plotZerosVsDepth <- function(object) {
    counts <- assay(object)
    metrics <- metrics(object)
    # using a logical matrix is much faster than a == 0L comparison
    present <- as(counts, "lgCMatrix")

    # The `counts == 0L` operation here blows up memory on large datasets
    df <- data.frame(
        dropout = (nrow(present) - Matrix::colSums(present)) / nrow(present),
        depth = Matrix::colSums(counts),
        fileName = metrics[["fileName"]])
    p <- ggplot(df,
           aes_(x = ~depth,
                y = ~dropout * 100L)) +
        geom_point(size = 0.8, alpha = 0.3) +
        scale_x_log10() +
        labs(x = "library size",
             y = "% genes zero")
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]])) {
        p <- p + facet_wrap(~fileName)
    }
    p
}



# Methods ====
#' @rdname plotZerosVsDepth
#' @export
setMethod("plotZerosVsDepth", "bcbioSingleCellANY", .plotZerosVsDepth)
