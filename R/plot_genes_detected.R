#' Genes detected plots
#'
#' @rdname plot_genes_detected
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param metrics Barcode metrics data frame
#' @param min_genes Recommended minimum gene count cutoff
#' @param max_genes Recommended maximum gene count cutoff
#'
#' @return ggplot2 object



#' @rdname plot_genes_detected
#' @description Boxplot
#' @export
plot_genes_detected_boxplot <- function(
    metrics,
    min_genes = get("min_genes", envir = parent.frame()),
    max_genes = get("max_genes", envir = parent.frame())) {
    boxplot <- metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~genes_detected,
                 fill = ~sample_name)
        ) +
        labs(title = "genes detected boxplot",
             x = "sample name",
             y = "genes per cell") +
        geom_boxplot() +
        geom_hline(color = warn_color,
                   yintercept = min_genes) +
        geom_hline(color = warn_color,
                   yintercept = max_genes) +
        geom_label(
            data = aggregate(genes_detected ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(genes_detected)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")
        ) +
        expand_limits(y = 0) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(boxplot)
}



#' @rdname plot_genes_detected
#' @description Histogram
#' @export
plot_genes_detected_histogram <- function(
    metrics,
    min_genes = get("min_genes", envir = parent.frame()),
    max_genes = get("max_genes", envir = parent.frame())) {
    histogram <- metrics %>%
        ggplot(
            aes_(x = ~genes_detected,
                 fill = ~sample_name)
        ) +
        labs(title = "genes detected histogram",
             x = "genes per cell") +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color,
                   xintercept = min_genes) +
        geom_vline(color = warn_color,
                   xintercept = max_genes) +
        expand_limits(x = 0) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(histogram)
}



#' @rdname plot_genes_detected
#' @description Show both plots
#'
#' @param ... Passthrough parameters
#'
#' @export
plot_genes_detected <- function(...) {
    show(plot_genes_detected_histogram(...))
    show(plot_genes_detected_boxplot(...))
}
