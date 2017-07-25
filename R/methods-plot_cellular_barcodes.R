#' Plot Cellular Barcode Distributions per Sample
#'
#' @rdname plot_cellular_barcodes
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @details A violin plot is a comact display of a continuous distribution. It
#'   is a blend of [geom_boxplot()] and [geom_density()]: a violin plot is a
#'   mirrored density plot displayed in the same way as a boxplot.
#'
#' @return [ggplot].



#' @rdname cellular_barcodes
#' @usage NULL
.cb_cutoff_line <- function(object) {
    metadata(object)[["cb_cutoff"]] %>%
        as.numeric %>%
        log10
}



#' @rdname cellular_barcodes
#' @usage NULL
.plot_cb_tbl <- function(object) {
    cellular_barcodes <- bcbio(object, "cellular_barcodes")
    if (is.null(cellular_barcodes)) {
        stop("Cellular barcode reads not saved in object")
    }
    interesting_group <- interesting_groups(object)[[1L]]
    sample_metadata <- sample_metadata(object) %>%
        tidy_select(unique(c("file_name",
                             "sample_id",
                             "sample_name",
                             interesting_group)))
    cellular_barcodes %>%
        .bind_cellular_barcodes %>%
        # Separate the cellular barcode with unique character (dot)
        mutate(cellular_barcode = str_replace(
            .data[["cellular_barcode"]],
            "_([acgt]{8}_[acgt]{8})$", ".\\1")) %>%
        separate_(col = "cellular_barcode",
                  into = c("sample_id", "cellular_barcode"),
                  sep = "\\.") %>%
        mutate(log10_reads = log10(.data[["reads"]]),
               reads = NULL) %>%
        # Only plot barcodes with at least 100 reads (log10 = 2)
        filter(.data[["log10_reads"]] > 2L) %>%
        left_join(sample_metadata, by = "sample_id")
}



#' @rdname plot_cellular_barcodes
#' @usage NULL
.plot_cb_raw_violin <- function(plot_cb_tbl, cb_cutoff_line) {
    ggplot(
        plot_cb_tbl,
        aes_(x = ~sample_name,
             y = ~log10_reads,
             fill = ~sample_name)) +
        facet_wrap(~file_name) +
        geom_violin(scale = "width") +
        geom_hline(alpha = qc_line_alpha,
                   color = qc_fail_color,
                   size = qc_line_size,
                   yintercept = 2L) +
        geom_hline(alpha = qc_line_alpha,
                   color = qc_pass_color,
                   size = qc_line_size,
                   yintercept = cb_cutoff_line) +
        labs(title = "raw violin",
             x = "",
             y = "log10 reads per cell") +
        coord_flip() +
        theme(legend.position = "none")
}



#' @rdname plot_cellular_barcodes
#' @usage NULL
.plot_cb_raw_histogram <- function(plot_cb_tbl, cb_cutoff_line) {
    ggplot(
        plot_cb_tbl,
        aes_(x = ~log10_reads,
             fill = ~sample_name)) +
        labs(title = "raw histogram",
             x = "log10 reads per cell",
             y = "") +
        facet_wrap(~file_name) +
        geom_histogram(bins = 100L) +
        scale_y_sqrt() +
        geom_vline(alpha = qc_line_alpha,
                   color = qc_fail_color,
                   size = qc_line_size,
                   xintercept = 2L) +
        geom_vline(alpha = qc_line_alpha,
                   color = qc_pass_color,
                   size = qc_line_size,
                   xintercept = cb_cutoff_line) +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              legend.position = "none")
}



#' @rdname plot_cellular_barcodes
#' @usage NULL
.plot_cb_proportional_histogram <- function(object) {
    cb_cutoff_line <- .cb_cutoff_line(object)
    ggplot(.proportional_cb(object),
           aes_(x = ~log10_reads_per_cell,
                y = ~proportion_of_cells * 100L,
                color = ~sample_name)) +
        facet_wrap(~file_name) +
        geom_line() +
        geom_vline(alpha = qc_line_alpha,
                   color = qc_fail_color,
                   size = qc_line_size,
                   xintercept = 2L) +
        geom_vline(alpha = qc_line_alpha,
                   color = qc_pass_color,
                   size = qc_line_size,
                   xintercept = cb_cutoff_line) +
        labs(title = "proportional histogram",
             x = "log10 reads per cell",
             y = "% of cells") +
        theme(legend.position = "bottom")
}



#' @rdname plot_cellular_barcodes
#' @export
setMethod("plot_cellular_barcodes", "bcbioSCDataSet", function(object) {
    # Use defined plot_cb_tbl here for improved speed
    plot_cb_tbl <- .plot_cb_tbl(object)
    cb_cutoff_line <- .cb_cutoff_line(object)
    ggdraw() +
        # Coordinates are relative to lower left corner
        draw_plot(
            .plot_cb_raw_violin(plot_cb_tbl, cb_cutoff_line),
            x = 0L, y = 0.7, width = 0.5, height = 0.3) +
        draw_plot(
            .plot_cb_raw_histogram(plot_cb_tbl, cb_cutoff_line),
            x = 0.5, y = 0.7, width = 0.5, height = 0.3) +
        draw_plot(
            .plot_cb_proportional_histogram(object),
            x = 0L, y = 0L, width = 1L, height = 0.7)
})
