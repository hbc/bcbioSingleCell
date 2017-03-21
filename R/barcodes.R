##' Load a barcode histogram file generated from bcbio-nextgen
##'
##' @param filename path to a barcode histogram file
##' @return dataframe of reads per barcode
##' @importFrom readr read_tsv
##' @author Rory Kirchner
##' @export
read_barcode_file <- function(filename) {
    readr::read_tsv(filename,
                    col_names = c("barcode", "count"),
                    progress = FALSE) %>%
        return
}

##' Plot a barcode histogram for a given sample
##'
##' @param filename path to a barcode histogram file
##' @param sample title for plot
##' @import ggplot2
##' @import scales
##' @importFrom graphics hist
##' @return a ggplot2 plot object
##' @author Rory Kirchner
##' @export
barcode_plot <- function(filename, sample) {
    bcs <- read_barcode_file(filename)
    bcs_hist <- graphics::hist(log10(bcs$count),
                               plot = FALSE, n = 50)

    fLog <- bcs_hist$count
    xLog <- bcs_hist$mids

    y <- fLog * (10^xLog) / sum(fLog * (10^xLog))

    ggplot2::qplot(10^xLog, y) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::ggtitle(sample) +
        ggplot2::scale_x_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(~10^.x))
        ) +
        ggplot2::xlab("number of reads assigned to a cell") +
        ggplot2::ylab("proportion of cells") %>%
        print
}
