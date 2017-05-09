#' Total counts plots
#'
#' @rdname plot_total_counts
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param metrics Barcode metrics data frame
#'
#' @return ggplot2 object



#' @rdname plot_total_counts
#' @description Histogram
#' @export
plot_total_counts_histogram <- function(metrics) {
    histogram <- metrics %>%
        ggplot(
            aes_(x = ~total_counts,
                 fill = ~sample_name)
        ) +
        labs(title = "total RNA read counts histogram",
             x = "counts per cell") +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        expand_limits(x = 0) +
        scale_y_sqrt() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(histogram)
}



#' @rdname plot_total_counts
#' @description Boxplot
#' @export
plot_total_counts_boxplot <- function(metrics) {
    boxplot <- metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~total_counts,
                 fill = ~sample_name)
        ) +
        labs(title = "total RNA read counts boxplot",
             x = "sample name",
             y = "counts per cell") +
        geom_boxplot() +
        geom_label(
            data = aggregate(total_counts ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(total_counts)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")
        ) +
        scale_y_log10() +
        expand_limits(y = 1) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(boxplot)
}



#' @rdname plot_total_counts
#' @description Show both plots
#'
#' @param ... Passthrough parameters
#'
#' @export
plot_total_counts <- function(...) {
    show(plot_total_counts_histogram(...))
    show(plot_total_counts_boxplot(...))
}
