##' Load a barcode histogram file generated from bcbio-nextgen
##'
##' @param filename path to a barcode histogram file
##' @return dataframe of reads per barcode
##' @importFrom readr read_tsv
##' @author Rory Kirchner
##' @export
read_barcode_file = function(filename) {
  return(read_tsv(filename, col_names = c("barcode", "count"), progress = FALSE))
}

##' Plot a barcode histogram for a given sample
##'
##' @param filename path to a barcode histogram file
##' @param sample title for plot
##' @import ggplot2
##' @import scales
##' @return a ggplot2 plot object
##' @author Rory Kirchner
##' @export
barcode_plot = function(filename, sample) {
  bcs = read_barcode_file(filename)
  bcs_hist = hist(log10(bcs$count), plot = FALSE, n = 50)
  fLog = bcs_hist$count
  xLog = bcs_hist$mids
  y = fLog * (10^xLog)/sum(fLog * (10^xLog))
  p = qplot(10^xLog, y) +
    geom_point() +
    theme_bw() +
    ggtitle(sample) +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    xlab("number of reads assigned to a cell") +
    ylab("proportion of cells")
  return(p)
}
