#' Total counts vs. genes detected plot
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#'
#' @param metrics Barcode metrics data frame
#' @param colorby column to color the points by
#' @return ggplot2 object
#' @export
plot_total_vs_detected <- function(metrics, colorby="sample") {
    plot <- metrics %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~total_counts,
                          y = ~genes_detected,
                          color = as.name(colorby))) +
        ggplot2::facet_wrap(~sample) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(method = "loess") +
        ggplot2::ggtitle("total counts vs. genes detected") +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90,
                                                hjust = 1)
        ) +
        ggplot2::xlab("counts per cell") +
        ggplot2::ylab("genes per cell")
    return(plot)
}
