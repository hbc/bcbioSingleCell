#' Plot Percentage of Zeros vs. Library Depth
#'
#' @rdname plotZerosVsDepth
#' @name plotZerosVsDepth
#'
#' @return [ggplot].
NULL



# Constructors ====
.plotZerosVsDepth <- function(object) {
    counts <- assay(object)
    zeros <- data.frame(
        dropout = Matrix::colSums(counts == 0L) / nrow(counts),
        depth = Matrix::colSums(counts))
    p <- ggplot(zeros,
           aes_(x = ~depth,
                y = ~dropout * 100L)) +
        geom_point(size = 0.8) +
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
setMethod("plotZerosVsDepth", "bcbioSCDataSet", .plotZerosVsDepth)



#' @rdname plotZerosVsDepth
#' @export
setMethod("plotZerosVsDepth", "bcbioSCSubset", .plotZerosVsDepth)
