#' Plot Percentage of Zeros vs. Library Depth
#'
#' This function helps us visualize the dropout rate.
#'
#' @name plotZerosVsDepth
#' @family Quality Control Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#' @param metrics Metrics `data.frame`.
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotZerosVsDepth(bcb_small)
#'
#' # seurat ====
#' plotZerosVsDepth(pbmc_small)
#' plotZerosVsDepth(seurat_small)
NULL



# Constructors =================================================================
#' @importFrom ggplot2 aes_string facet_wrap geom_point ggplot labs
#'   scale_x_log10
.plotZerosVsDepth <- function(object) {
    counts <- counts(object)
    metrics <- metrics(object)

    # Using a logical matrix is faster and more memory efficient
    present <- counts %>%
        # Ensure dgTMatrix gets coereced (e.g. pbmc_small)
        as("dgCMatrix") %>%
        as("lgCMatrix")

    # Add dropout rate and depth
    metrics <- mutate(
        metrics,
        "dropout" = (nrow(present) - Matrix::colSums(present)) / nrow(present),
        "depth" = Matrix::colSums(counts)
    )

    p <- ggplot(
        data = metrics,
        mapping = aes_string(
            x = "depth",
            y = "dropout",
            color = "sampleName"
        )
    ) +
        geom_point(size = 0.8, alpha = 0.8) +
        scale_x_log10() +
        labs(x = "library size (depth)", y = "dropout rate")

    # Wrap aggregated samples
    facets <- NULL
    if (isTRUE(.checkAggregate(metrics))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    p
}



# Methods ======================================================================
#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("bcbioSingleCell"),
    .plotZerosVsDepth
)



#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("seurat"),
    .plotZerosVsDepth
)
