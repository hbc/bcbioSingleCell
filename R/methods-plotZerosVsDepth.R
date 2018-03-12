#' Plot Percentage of Zeros vs. Library Depth
#'
#' This function helps us visualize the dropout rate.
#'
#' @name plotZerosVsDepth
#' @family Quality Control Metrics
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param metrics Metrics `data.frame`.
#'
#' @return `ggplot`.
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' plotZerosVsDepth(bcb)
#'
#' # dgCMatrix ====
#' counts <- counts(bcb)
#' metrics <- metrics(bcb)
#' plotZerosVsDepth(counts, metrics = metrics)
#'
#' # seurat ====
#' plotZerosVsDepth(pbmc_small)
#' plotZerosVsDepth(seurat)
NULL



# Constructors =================================================================
#' @importFrom ggplot2 aes_string facet_wrap geom_point ggplot labs
#'   scale_x_log10
.plotZerosVsDepth <- function(object, metrics) {
    # Using a logical matrix is fast and memory efficient
    present <- as(object, "lgCMatrix")
    df <- data.frame(
        dropout = (nrow(present) - Matrix::colSums(present)) / nrow(present),
        depth = Matrix::colSums(object),
        description = metrics[["description"]]
    )
    p <- ggplot(
        df,
        mapping = aes_string(x = "depth", y = "dropout")
    ) +
        geom_point(size = 0.8, alpha = 0.3) +
        scale_x_log10() +
        labs(x = "library size (depth)", y = "dropout rate")

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}



# Methods ======================================================================
#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("bcbioSingleCell"),
    function(object) {
        .plotZerosVsDepth(
            object = counts(object),
            metrics = metrics(object)
        )
    }
)



#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("dgCMatrix"),
    .plotZerosVsDepth
)



#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("matrix"),
    .plotZerosVsDepth
)



#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("seurat"),
    function(object) {
        .plotZerosVsDepth(
            object = counts(object),
            metrics = metrics(object)
        )
    }
)
