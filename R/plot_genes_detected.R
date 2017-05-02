#' Genes detected plots
#'
#' @rdname plot_genes_detected
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param metrics Barcode metrics data frame
#'
#' @return ggplot2 object



#' @rdname plot_genes_detected
#' @description Boxplot
#' @export
plot_genes_detected_boxplot <- function(metrics) {
    boxplot <- metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~genes_detected,
                 fill = ~sample_name)
        ) +
        geom_boxplot() +
        geom_label(
            data = aggregate(genes_detected ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(genes_detected))
        ) +
        ggtitle("genes detected boxplot") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        ) +
        ylab("genes per cell")
    return(boxplot)
}



#' @rdname plot_genes_detected
#' @description Histogram
#' @export
plot_genes_detected_histogram <- function(metrics) {
    histogram <- metrics %>%
        ggplot(
            aes_(x = ~genes_detected,
                 fill = ~sample_name)
        ) +
        facet_wrap(~sample_name) +
        geom_histogram() +
        ggtitle("genes detected histogram") +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none") +
        xlab("genes per cell")
    return(histogram)
}



#' @rdname plot_genes_detected
#' @description Show both plots
#' @export
plot_genes_detected <- function(metrics) {
    show(plot_genes_detected_histogram(metrics))
    show(plot_genes_detected_boxplot(metrics))
}
