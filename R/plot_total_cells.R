#' Total cells plot
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param metrics Barcode metrics data frame
#'
#' @return ggplot2 object
#' @export
plot_total_cells <- function(metrics) {
    plot <- metrics %>%
        group_by_(.dots = "sample_name") %>%
        summarize_(total_cells = ~n()) %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~total_cells,
                 fill = ~sample_name)
        ) +
        geom_bar(stat = "identity") +
        geom_text(vjust = -0.5,
                  aes_(label = ~total_cells)
        ) +
        ggtitle("total number of cells") +
        theme(legend.position = "none") +
        ylab("cell count (w/ barcode depth cutoff)")
    return(plot)
}
