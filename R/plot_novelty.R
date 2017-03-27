#' Novelty plot
#'
#' log detected per count
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @import ggplot2
#' @importFrom stats aggregate median
#'
#' @param metrics Barcode metrics data frame
#'
#' @export
plot_novelty <- function(metrics) {
    histogram <- metrics %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~log_detected_per_count,
                          fill = ~sample)
        ) +
        ggplot2::facet_wrap(~sample) +
        ggplot2::geom_histogram() +
        ggplot2::ggtitle("novelty histogram") +
        ggplot2::scale_y_sqrt() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab("log genes detected per count")

    boxplot <- metrics %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~sample,
                          y = ~log_detected_per_count,
                          fill = ~sample)
        ) +
        ggplot2::geom_boxplot() +
        ggplot2::ggtitle("novelty boxplot") +
        ggplot2::geom_label(
            data = stats::aggregate(log_detected_per_count ~ sample,
                                    metrics,
                                    stats::median),
            aes_(label = ~round(log_detected_per_count, digits = 2))
        ) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ylab("log genes detected per count") +
        ggplot2::expand_limits(y = 0)

    print(histogram)
    print(boxplot)
}
