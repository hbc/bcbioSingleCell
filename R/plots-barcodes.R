#' Plot cellular barcode distributions per sample
#'
#' @details A violin plot is a compact display of a continuous distribution. It
#'   is a blend of `geom_boxplot()` and `geom_density()`: a violin plot is a
#'   mirrored density plot displayed in the same way as a boxplot.
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param object [bcbioSCDataSet].

#' @export



#' @rdname plot_cellular_barcodes
#' @usage NULL
.plot_cb_df <- function(object) {
    message("Reformatting cellular barcodes for ggplot")
    cellular_barcodes <- bcbio(object, "cellular_barcodes")
    if (is.null(cellular_barcodes)) {
        return(NULL)
    }
    interesting_group <- interesting_groups(object)[[1L]]
    sample_metadata <- sample_metadata(object) %>%
        tidy_select(unique(c("file_name",
                             "sample_id",
                             "sample_name",
                             interesting_group)))
    cellular_barcodes %>%
        .bind_cellular_barcodes %>%
        separate_(col = "cellular_barcode",
                  into = c("sample_id", "cellular_barcode"),
                  sep = ":") %>%
        mutate(log10_reads = log10(.data[["reads"]]),
               reads = NULL) %>%
        filter(.data[["log10_reads"]] > 0L) %>%
        left_join(sample_metadata, by = "sample_id")
}



#' @rdname plot_cellular_barcodes
#' @usage NULL
.plot_cb_violin <- function(plot_cb_df) {
    ggplot(
        plot_cb_df,
        aes_(x = ~sample_name,
             y = ~log10_reads,
             fill = ~sample_name)) +
        facet_wrap(~file_name) +
        geom_violin(scale = "width") +
        labs(title = "cellular barcode violin plot",
             x = "sample name",
             y = "log10 reads per cell") +
        coord_flip()
}



#' @rdname plot_cellular_barcodes
#' @usage NULL
.plot_cb_histogram <- function(plot_cb_df) {
    ggplot(
        plot_cb_df,
        aes_(x = ~log10_reads,
             fill = ~sample_name)) +
        labs(title = "cellular barcode histogram",
             x = "log10 reads per cell") +
        facet_wrap(~file_name) +
        geom_histogram(bins = 200L) +
        scale_y_sqrt()
}



#' @rdname plot_cellular_barcodes
#' @usage NULL
.plot_cb_histogram2 <- function(object) {
    cellular_barcodes <- bcbio(object, "cellular_barcodes")
    list <- lapply(seq_along(cellular_barcodes), function(a) {
        bcs <- cellular_barcodes[[a]] %>%
            mutate(log10_reads = log10(.data[["reads"]]))
        bcs_hist <- hist(bcs[["log10_reads"]], n = 100L, plot = FALSE)
        f_log <- bcs_hist[["counts"]]
        x_log <-  bcs_hist[["mids"]]
        x <- x_log
        y <- f_log * (10L ^ x_log) / sum(f_log * (10L ^ x_log))
        data.frame(
            sample_id = names(cellular_barcodes)[[a]],
            x = x,
            y = y)
    }) %>%
        set_names(names(cellular_barcodes))
    bind_rows(list) %>%
        left_join(sample_metadata(object), by = "sample_id") %>%
        ggplot(
            aes_(x = ~x,
                 y = ~y,
                 color = ~sample_name)) +
        facet_wrap(~file_name) +
        geom_line() +
        labs(title = "cellular barcode histogram",
             x = "log10 reads per cell",
             y = "proportion of cells")
}



#' @rdname plot_cellular_barcodes
#' @export
plot_cellular_barcodes <- function(object) {
    plot_cb_df <- .plot_cb_df(object)
    if (is.null(plot_cb_df)) {
        return(NULL)
    }
    plot_grid(.plot_cb_violin(plot_cb_df),
              .plot_cb_histogram(plot_cb_df),
              .plot_cb_histogram2(object),
              labels = "auto",
              nrow = 3L)
}
