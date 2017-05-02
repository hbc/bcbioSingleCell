#' Mitochondrial count plots
#'
#' @rdname plot_mito_counts
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param metrics Barcode metrics data frame
#'
#' @return ggplot2 object



#' @rdname plot_mito_counts
#' @description Histogram
#' @export
plot_mito_counts_histogram <- function(metrics) {
    metrics <- mutate(metrics, percent_mito = .data$percent_mito * 100)
    histogram <- metrics %>%
        ggplot(
            aes_(x = ~percent_mito,
                 fill = ~sample_name)) +
        facet_wrap(~sample_name) +
        geom_histogram() +
        geom_vline(color = warn_color, xintercept = 20) +
        ggtitle("mitochondrial gene abundance histogram") +
        scale_y_sqrt() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        ) +
        xlab("% mitochondrial") +
        xlim(0, 100)
    return(histogram)
}



#' @rdname plot_mito_counts
#' @description Boxplot
#' @export
plot_mito_counts_boxplot <- function(metrics) {
    metrics <- mutate(metrics, percent_mito = .data$percent_mito * 100)
    boxplot <- metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~percent_mito,
                 fill = ~sample_name)
        ) +
        geom_boxplot() +
        geom_hline(color = warn_color, yintercept = 20) +
        ggtitle("mitochondrial gene abundance boxplot") +
        geom_label(
            data = aggregate(percent_mito ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(percent_mito, digits = 2))
        ) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        ) +
        ylab("% mitochondrial") +
        expand_limits(y = 0)
    return(boxplot)
}



#' @rdname plot_mito_counts
#' @description Scatterplot
#' @export
plot_mito_counts_scatterplot <- function(metrics) {
    scatterplot <- metrics %>%
        ggplot(
            aes_(x = ~coding_counts,
                 y = ~mito_counts,
                 color = ~sample_name)
        ) +
        facet_wrap(~sample_name) +
        ggtitle("mitochondrial gene abundance scatterplot") +
        geom_point() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        ) +
        xlab("counts in mitochondrial genes") +
        # xlim(0, 50000) +
        ylab("counts in coding genes")
    return(scatterplot)
}



#' @rdname plot_mito_counts
#' @description Show both plots
#' @export
plot_mito_counts <- function(metrics) {
    show(plot_mito_counts_histogram(metrics))
    show(plot_mito_counts_boxplot(metrics))
    show(plot_mito_counts_scatterplot(metrics))
}
