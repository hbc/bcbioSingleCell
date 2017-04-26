#' Novelty plots (log detected per count)
#'
#' @rdname plot_novelty
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param metrics Barcode metrics data frame
#'
#' @return ggplot2 object



#' @rdname plot_novelty
#' @description Histogram
#' @export
plot_novelty_histogram <- function(metrics) {
    histogram <- metrics %>%
        ggplot(
            aes_(x = ~log10_detected_per_count,
                 fill = ~sample_name)
        ) +
        facet_wrap(~sample_name) +
        geom_histogram() +
        ggtitle("novelty histogram") +
        scale_y_sqrt() +
        theme(legend.position = "none") +
        xlab("log genes detected per count")
    return(histogram)
}



#' @rdname plot_novelty
#' @description Boxplot
#' @export
plot_novelty_boxplot <- function(metrics) {
    boxplot <- metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~log10_detected_per_count,
                 fill = ~sample_name)
        ) +
        geom_boxplot() +
        ggtitle("novelty boxplot") +
        geom_label(
            data = aggregate(log10_detected_per_count ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(log10_detected_per_count, digits = 2))
        ) +
        theme(legend.position = "none") +
        ylab("log genes detected per count") +
        expand_limits(y = 0)
    return(boxplot)
}



#' @rdname plot_novelty
#' @description Show both plots
#' @export
plot_novelty <- function(metrics) {
    show(plot_novelty_histogram(metrics))
    show(plot_novelty_boxplot(metrics))
}
