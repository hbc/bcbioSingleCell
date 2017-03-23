#' Total counts plots
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats aggregate median
#'
#' @param metrics Barcode metrics data frame
#'
#' @export
plot_total_counts <- function(metrics) {
    histogram <- metrics %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~total_counts,
                          fill = ~sample)
        ) +
        ggplot2::facet_wrap(~sample) +
        ggplot2::geom_histogram() +
        ggplot2::ggtitle("total RNA read counts histogram") +
        ggplot2::scale_x_log10() +
        ggplot2::scale_y_log10() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
            legend.position = "none"
        ) +
        ggplot2::xlab("counts per cell")

    boxplot <- metrics %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~sample,
                          y = ~total_counts,
                          fill = ~sample)
        ) +
        ggplot2::geom_boxplot() +
        ggplot2::geom_label(
            data = stats::aggregate(total_counts ~ sample,
                                    metrics,
                                    stats::median),
            ggplot2::aes_(label = ~round(total_counts))
        ) +
        ggplot2::ggtitle("total RNA read counts boxplot") +
        ggplot2::scale_y_log10() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ylab("counts per cell")

    print(histogram)
    print(boxplot)
}
