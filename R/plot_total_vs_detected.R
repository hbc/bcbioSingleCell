#' Total counts vs. genes detected plot
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param metrics Barcode metrics data frame
#' @param colorby column to color the points by
#'
#' @return ggplot2 object
#' @export
plot_total_vs_detected <- function(metrics, colorby = "sample") {
    plot <- metrics %>%
        ggplot(
            aes_(x = ~total_counts,
                 y = ~genes_detected,
                 color = as.name(colorby))) +
        facet_wrap(~sample) +
        geom_point() +
        geom_smooth(method = "loess") +
        ggtitle("total counts vs. genes detected") +
        scale_x_log10() +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        xlab("counts per cell") +
        ylab("genes per cell")
    return(plot)
}
