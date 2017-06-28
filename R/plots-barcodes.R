# Cellular barcode distributions ====
# Histogram function updated to use stored barcode counts rather than by
# accessing the files directly on the cluster. This approach allows us to
# generate multiple types of plots more easily, such as the new violin plot
# method.
.plot_cb_histogram <- function(bcb) {
    barcodes <- bcbio(bcb, "cellular_barcodes")
    message("Generating barcode histograms")
    pblapply(seq_along(barcodes), function(a) {
        title <- paste0(
            names(barcodes)[[a]],
            " (",
            sample_metadata(bcb)[name, "sample_name"],
            ")")
        bcs <- barcodes[[a]] %>%
            mutate(log10_reads = log10(.data[["reads"]]))
        bcs_hist <- hist(bcs[["log10_reads"]], n = 50L, plot = FALSE)
        f_log <- bcs_hist[["counts"]]
        x_log <- bcs_hist[["mids"]]
        qplot(x = 10L ^ x_log,
              y = f_log * (10L ^ x_log) / sum(f_log * (10L ^ x_log))) +
            geom_point() +
            geom_line() +
            scale_x_log10(
                breaks = trans_breaks("log10", function(x) 10L ^ x),
                labels = trans_format("log10", math_format(~10L ^ .x))) +
            labs(title = title,
                 x = "reads per cell",
                 y = "proportion of cells")
    }) %>% invisible
}

.plot_cb_violin <- function(bcb) {
    barcodes <- bcbio(bcb, "cellular_barcodes") %>%
        .bind_cellular_barcodes %>%
        separate_(col = "rowname",
                  into = c("sample_id", "cellular_barcode"),
                  sep = ":") %>%
        mutate(cellular_barcode = NULL,
               log10_reads = log10(.data[["reads"]])) %>%
        left_join(sample_metadata(bcb), by = "sample_id")
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(barcodes,
           aes_(x = ~sample_name,
                y = ~log10_reads,
                fill = as.name(interesting_group))) +
        facet_wrap(~file_name) +
        geom_violin(scale = "width") +
        labs(title = "cellular barcode violin plot",
             x = "sample name",
             y = "log10 reads per cellular barcode") +
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
    show(.plot_cb_violin(bcb))
    show(.plot_cb_histogram(bcb))
}
