#' Total cells plot
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#' @keywords plot qc

#' @import dplyr
#' @import ggplot2
#'
#' @param metrics Barcode metrics data frame
#'
#' @export
total_cells_plot <- function(metrics) {
    plot <- metrics %>%
        dplyr::group_by_(.dots = "sample") %>%
        dplyr::summarize_(total_cells = ~n()) %>%
        ggplot2::ggplot(
            ggplot2::aes_(x = ~sample,
                          y = ~total_cells,
                          fill = ~sample)
        ) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_text(
            ggplot2::aes_(label = ~total_cells)
        ) +
        ggplot2::ggtitle("total number of cells") +
        ggplot2::theme(legend.position = "none") +
        ggplot2::ylab("cell count (w/ barcode cutoff)")
    print(plot)
}
