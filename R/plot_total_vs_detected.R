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
plot_total_vs_detected <- function(metrics, colorby = "sample_name") {
    plot <- metrics %>%
        ggplot(
            aes_(x = ~total_counts,
                 y = ~genes_detected,
                 color = as.name(colorby))
        ) +
        labs(title = "total counts vs. genes detected",
             x = "counts per cell",
             y = "genes per cell") +
        facet_wrap(~sample_name) +
        geom_point() +
        geom_smooth(method = "loess") +
        scale_x_log10() +
        scale_y_log10() +
        # expand_limits(x = 1, y = 1) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(plot)
}
