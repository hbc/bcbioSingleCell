# Cellular barcode distributions ====
# Histogram function updated to use stored barcode counts rather than by
# accessing the files directly on the cluster. This approach allows us to
# generate multiple types of plots more easily, such as the new violin plot
# method.
.plot_cb_histogram <- function(bcb) {
    barcodes <- bcb$barcodes
    message("Generating barcode histograms...")
    pbmclapply(seq_along(barcodes), function(a) {
        name <- names(bcb$barcodes[a])
        sample_name <- bcb$metadata[name, "sample_name"]
        title <- paste(sample_name, name, sep = " : ")
        bcs <- bcb$barcodes[a] %>%
            as.data.frame %>%
            rownames_to_column %>%
            set_colnames(c("cellular_barcode", "reads")) %>%
            mutate(log10_reads = log10(.data$reads))
        bcs_hist <- hist(bcs$log10_reads, plot = FALSE, n = 50)
        fLog <- bcs_hist$count
        xLog <- bcs_hist$mids
        y <- fLog * (10^xLog) / sum(fLog * (10^xLog))
        qplot(10^xLog, y) +
            geom_point() +
            geom_line() +
            ggtitle(title) +
            scale_x_log10(
                breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(~10^.x))) +
            xlab("number of reads assigned to a cell") +
            ylab("proportion of cells")
    }) %>% invisible
}

.plot_cb_violin <- function(bcb) {
    title <- deparse(substitute(bcb))
    barcodes <- unlist_barcodes(bcb)
    message("Generating violin plot...")
    ggplot(barcodes,
           aes_(x = ~sample_name,
                y = ~log10_reads,
                fill = ~sample_name)) +
        geom_violin(scale = "width") +
        labs(title = paste(title, "violin plot"),
             x = "sample name",
             y = "log10 reads per cellular barcode") +
        theme(legend.position = "none") +
        coord_flip()
}

#' Plot cellular barcode distributions per sample
#'
#' @details A violin plot is a compact display of a continuous distribution. It
#'   is a blend of `geom_boxplot()` and `geom_density()`: a violin plot is a
#'   mirrored density plot displayed in the same way as a boxplot.
#'
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param bcb [bcbioSCDataSet].

#' @export
plot_barcodes <- function(bcb) {
    .plot_cb_violin(bcb) %>% show
    .plot_cb_histogram(bcb) %>% show
}
