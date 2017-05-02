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
#' @param file_name Barcode histogram file
#' @param sample_name Sample name (title for plot)
#'
#' @export
plot_barcode <- function(file_name, sample_name = NULL) {
    # Get the sample name from the file name by default
    if (is.null(sample_name)) {
        sample_name <- gsub("-barcodes\\.tsv$", "", basename(file_name))
    }

    bcs <- read_barcode_file(file_name)
    bcs_hist <- hist(log10(bcs$count), plot = FALSE, n = 50)

    fLog <- bcs_hist$count
    xLog <- bcs_hist$mids

    y <- fLog * (10^xLog) / sum(fLog * (10^xLog))

    plot <- qplot(10^xLog, y) +
        geom_point() +
        geom_line() +
        ggtitle(sample_name) +
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
    sample_barcodes <- names(run$sample_dirs)
    file_names <- file.path(
        run$sample_dirs, paste0(sample_barcodes, "-barcodes.tsv"))
    if (!all(file.exists(file_names))) {
        stop("Could not locate barcode TSV files")
    }
    metadata <- run$metadata[sample_barcodes,
                             c("sample_barcode", "sample_name")]
    names(file_names) <- paste(metadata$sample_name,
                               names(file_names),
                               sep = " : ")
    # Order the file names by sample name
    file_names <- file_names[order(names(file_names))]
    # Iterate over the files and show barcode plots
    lapply(seq_along(file_names), function(a) {
        show(plot_barcode(file_names[a], names(file_names)[a]))
    }) %>% invisible
}
