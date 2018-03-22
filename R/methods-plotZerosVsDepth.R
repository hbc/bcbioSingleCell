#' Plot Percentage of Zeros vs. Library Depth
#'
#' This function helps us visualize the dropout rate.
#'
#' @name plotZerosVsDepth
#' @family Quality Control Metrics
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
#' # dgCMatrix ====
#' counts <- counts(bcb_small)
#' metrics <- metrics(bcb_small)
#' plotZerosVsDepth(counts, metrics = metrics)
#'
#' # seurat ====
#' plotZerosVsDepth(pbmc_small)
#' plotZerosVsDepth(seurat_small)
NULL



# Constructors =================================================================
.plotZerosVsDepth <- function(object, metrics) {
    # Using a logical matrix is faster and more memory efficient
    present <- object %>%
        # Ensure dgTMatrix gets coereced (e.g. pbmc_small)
        as("dgCMatrix") %>%
        as("lgCMatrix")
    data <- tibble(
        "dropout" = (nrow(present) - Matrix::colSums(present)) / nrow(present),
        "depth" = Matrix::colSums(object),
        "sampleName" = metrics[["sampleID"]]
    )
    ggplot(
        data = data,
        mapping = aes_string(
            x = "depth",
            y = "dropout",
            color = "sampleName"
        )
    ) +
        geom_point(size = 0.8, alpha = 0.8) +
        scale_x_log10() +
        labs(x = "library size (depth)", y = "dropout rate")
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
