#' Total counts vs. genes detected plot
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#' @keywords plot qc
#'
#' @import dplyr
#' @import ggplot2
#' @importFrom stats setNames
#'
#' @param metrics Barcode metrics data frame
#'
#' @export
mito_counts_plot <- function(metrics) {
    histogram <- metrics %>%
        dplyr::mutate_(.dots = stats::setNames(
            list(~percent_mito * 100),
            "percent_mito"
        )) %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~percent_mito,
                          fill = ~sample)) +
        ggplot2::facet_wrap(~sample) +
        ggplot2::geom_histogram() +
        ggplot2::ggtitle("mitochondrial gene abundance histogram") +
        ggplot2::scale_y_sqrt() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::xlab("% mitochondrial")

    scatterplot <- metrics %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~coding_counts,
                          y = ~mito_counts,
                          color = ~sample)
        ) +
        ggplot2::facet_wrap(~sample) +
        ggplot2::ggtitle("mitochondrial gene abundance scatterplot") +
        ggplot2::geom_point() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90,
                                                hjust = 1),
            legend.position = "none"
        ) +
        ggplot2::xlab("counts in mitochondrial genes") +
        ggplot2::xlim(0, 50000) +
        ggplot2::ylab("counts in coding genes")

    print(histogram)
    print(scatterplot)
}


#' @rdname mito_counts_plot
#' @export
mito_counts_plots <- mito_counts_plot
