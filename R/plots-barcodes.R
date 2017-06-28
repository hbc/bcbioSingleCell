.plot_cellular_barcode_histogram <- function(bcb) {
    cellular_barcodes <- bcbio(bcb, "cellular_barcodes")
    if (is.null(cellular_barcodes)) {
        return(NULL)
    }
    message("Calculating barcode histograms")
    list <- pblapply(seq_along(cellular_barcodes), function(a) {
        bcs <- cellular_barcodes[[a]] %>%
            mutate(log10_reads = log10(.data[["reads"]]))
        bcs_hist <- hist(bcs[["log10_reads"]], n = 100L, plot = FALSE)
        f_log <- bcs_hist[["counts"]]
        x_log <-  bcs_hist[["mids"]]
        x <- 10L ^ x_log
        y <- f_log * (10L ^ x_log) / sum(f_log * (10L ^ x_log))
        data.frame(
            sample_id = names(cellular_barcodes)[[a]],
            x = x,
            y = y)
    }
    ) %>%
        set_names(names(barcodes))
    cb_bind <- bind_rows(list) %>%
        left_join(sample_metadata(bcb), by = "sample_id")
    ggplot(
        cb_bind,
        aes_(x = ~x,
             y = ~y,
             color = ~sample_name)) +
        facet_wrap(~file_name) +
        geom_point(size = 0.5) +
        geom_line() +
        scale_x_log10() +
        labs(title = "cellular barcode histogram",
             x = "log10 reads per cell",
             y = "proportion of cells")
}



.plot_cellular_barcode_violin <- function(bcb) {
    cellular_barcodes <- bcbio(bcb, "cellular_barcodes")
    if (is.null(cellular_barcodes)) {
        return(NULL)
    }
    cellular_barcodes <- cellular_barcodes %>%
        .bind_cellular_barcodes %>%
        separate_(col = "rowname",
                  into = c("sample_id", "cellular_barcode"),
                  sep = ":") %>%
        mutate(cellular_barcode = NULL,
               log10_reads = log10(.data[["reads"]])) %>%
        left_join(sample_metadata(bcb), by = "sample_id")
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(cellular_barcodes,
           aes_(x = ~sample_name,
                y = ~log10_reads,
                fill = as.name(interesting_group))) +
        facet_wrap(~file_name) +
        geom_violin(scale = "width") +
        labs(title = "cellular barcode violin plot",
             x = "sample name",
             y = "log10 reads per cell") +
        coord_flip()
}



#' Plot cellular barcode distributions per sample
#'
#' @details A violin plot is a compact display of a continuous distribution. It
#'   is a blend of `geom_boxplot()` and `geom_density()`: a violin plot is a
#'   mirrored density plot displayed in the same way as a boxplot.
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param bcb [bcbioSCDataSet].

#' @export
plot_cellular_barcodes <- function(bcb) {
    show(.plot_cellular_barcode_histogram(bcb))
    show(.plot_cellular_barcode_violin(bcb))
}
