##' Plot a barcode histogram for a given sample
##'
##' @param filename path to a barcode histogram file
##' @param sample title for plot
##' @import ggplot2
##' @importFrom graphics hist
##' @importFrom scales math_format trans_breaks trans_format
##' @keywords internal
##' @author Rory Kirchner
##' @author Michael Steinbaugh
##' @export
plot_barcode <- function(filename, sample = NULL) {
    # Get the sample name from the filename by default
    if (is.null(sample)) {
        sample <- gsub("-barcodes\\.tsv$", "", basename(filename))
    }

    bcs <- read_barcode_file(filename)
    bcs_hist <- graphics::hist(log10(bcs$count),
                               plot = FALSE, n = 50)

    fLog <- bcs_hist$count
    xLog <- bcs_hist$mids

    y <- fLog * (10^xLog) / sum(fLog * (10^xLog))

    p <- ggplot2::qplot(10^xLog, y) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::ggtitle(sample) +
        ggplot2::scale_x_log10(
            breaks = scales::trans_breaks("log10", function(x) 10^x),
            labels = scales::trans_format("log10", scales::math_format(~10^.x))
        ) +
        ggplot2::xlab("number of reads assigned to a cell") +
        ggplot2::ylab("proportion of cells")

    print(p)
}
