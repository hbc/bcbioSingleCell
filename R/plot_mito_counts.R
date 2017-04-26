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
    histogram <- metrics %>%
        mutate_(.dots = set_names(
            list(~percent_mito * 100),
            "percent_mito"
        )) %>%
        ggplot(
            aes_(x = ~percent_mito,
                 fill = ~sample_name)) +
        facet_wrap(~sample_name) +
        geom_histogram() +
        ggtitle("mitochondrial gene abundance histogram") +
        scale_y_sqrt() +
        theme(legend.position = "none") +
        xlab("% mitochondrial")
    return(histogram)
}



#' @rdname plot_mito_counts
#' @description Boxplot
#' @export
plot_mito_counts_boxplot <- function(metrics) {
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
            axis.text.x = element_text(angle = 90,
                                       hjust = 1),
            legend.position = "none"
        ) +
        xlab("counts in mitochondrial genes") +
        xlim(0, 50000) +
        ylab("counts in coding genes")
    return(scatterplot)
}



#' @rdname plot_mito_counts
#' @description Show both plots
#' @export
plot_mito_counts <- function(metrics) {
    show(plot_mito_counts_histogram(metrics))
    show(plot_mito_counts_boxplot(metrics))
}
