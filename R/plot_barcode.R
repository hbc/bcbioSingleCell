#' Plot sample barcodes
#'
#' @rdname plot_barcode
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh



#' @rdname plot_barcode
#' @description Plot a barcode histogram for a given sample
#' @keywords internal
#'
#' @param file_name path to a barcode histogram file
#' @param sample title for plot
#'
#' @export
plot_barcode <- function(file_name, sample = NULL) {
    # Get the sample name from the file name by default
    if (is.null(sample)) {
        sample <- gsub("-barcodes\\.tsv$", "", basename(file_name))
    }

    bcs <- read_barcode_file(file_name)
    bcs_hist <- hist(log10(bcs$count), plot = FALSE, n = 50)

    fLog <- bcs_hist$count
    xLog <- bcs_hist$mids

    y <- fLog * (10^xLog) / sum(fLog * (10^xLog))

    plot <- qplot(10^xLog, y) +
        geom_point() +
        geom_line() +
        ggtitle(sample) +
        scale_x_log10(
            breaks = trans_breaks("log10", function(x) 10^x),
            labels = trans_format("log10", math_format(~10^.x))
        ) +
        xlab("number of reads assigned to a cell") +
        ylab("proportion of cells")

    return(plot)
}



#' @rdname plot_barcode
#' @description Plot all single cell sample barcodes
#'
#' @param run \code{bcbio-nextgen} run
#'
#' @export
plot_barcodes <- function(run) {
    files <- list.files(
        run$final_dir,
        pattern = "*-barcodes.tsv",
        recursive = TRUE,
        include.dirs = TRUE,
        full.names = TRUE)
    lapply(seq_along(files), function(a) {
        show(plot_barcode(files[a]))
    }) %>% invisible
}
