#' Genes detected boxplot
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats aggregate median
#'
#' @param metrics Barcode metrics data frame
#' @return ggplot2 object
#' @export
plot_genes_detected_boxplot <- function(metrics) {
    boxplot <- metrics %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~sample,
                          y = ~genes_detected,
                          fill = ~sample)
        ) +
        ggplot2::geom_boxplot() +
        ggplot2::geom_label(
            data = stats::aggregate(genes_detected ~ sample,
                                    metrics,
                                    stats::median),
            ggplot2::aes_(label = ~round(genes_detected))
        ) +
        ggtitle("genes detected boxplot") +
        theme(legend.position = "none") +
        ylab("genes per cell")
    return(boxplot)
}

#' Genes detected histogram
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats aggregate median
#' @param metrics
#' @return ggplot2 object
plot_genes_detected_histogram <- function(metrics) {
    histogram <- metrics %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~genes_detected,
                          fill = ~sample)
        ) +
        ggplot2::facet_wrap(~sample) +
        ggplot2::geom_histogram() +
        ggplot2::ggtitle("genes detected histogram") +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90,
                                                hjust = 1),
            legend.position = "none") +
        ggplot2::xlab("genes per cell")
    return(histogram)
}
