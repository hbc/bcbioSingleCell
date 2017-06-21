#' Generate barcode metrics summary
#'
#' @rdname metrics
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run [bcbioSCDataSet].
#'
#' @return Tibble grouped by sample name.
.metrics <- function(run) {
    # Check for [Matrix::colSums()] methods support
    if (!"colSums,dgCMatrix-method" %in% methods(colSums)) {
        stop("dgCMatrix not supported in `colSums()`")
    }

    counts <- run$counts
    metadata <- run$metadata

    # Ensembl annotations, with broad class definitions (coding, mito)
    ensembl <- run$ensembl %>%
        ungroup %>%
        mutate(ensembl_transcript_id = NULL) %>%
        distinct
    coding <- ensembl %>%
        filter(.data$broad_class == "coding") %>%
        tidy_select(.data$external_gene_name) %>%
        .[[1]] %>% unique %>% sort
    mito <- ensembl %>%
        filter(.data$broad_class == "mito") %>%
        tidy_select(.data$external_gene_name) %>%
        .[[1]] %>% unique %>% sort

    data.frame(
        # Unique: `file_name` + `sample_barcode` + `cellular_barcode`
        unique = colnames(counts),
        total_counts = Matrix::colSums(counts),
        genes_detected = Matrix::colSums(counts > 0),
        coding_counts = Matrix::colSums(
            counts[rownames(counts) %in% coding, ]),
        mito_counts = Matrix::colSums(
            counts[rownames(counts) %in% mito, ])) %>%
        # Separate the barcodes, later used to join metadata
        separate_("unique",
                  c("sample_barcode", "cellular_barcode"),
                  sep = ":",
                  remove = FALSE) %>%
        # Summary statistics
        mutate(unique = NULL,
               log10_detected_per_count = log10(.data$genes_detected) /
                   log10(.data$total_counts),
               mito_ratio = .data$mito_counts / .data$total_counts) %>%
        left_join(metadata[, c("sample_barcode", "sample_name")],
                  by = "sample_barcode") %>%
        # Select sample name first
        tidy_select(.data$sample_name,
                    .data$sample_barcode,
                    .data$cellular_barcode,
                    everything()) %>%
        group_by(!!!syms(c("sample_name",  "sample_barcode"))) %>%
        arrange(desc(!!sym("total_counts")), .by_group = TRUE) %>%
        ungroup %>%
        mutate(rowname = paste(.data$sample_barcode,
                               .data$cellular_barcode,
                               sep = ":")) %>%
        as.data.frame %>%
        column_to_rownames %>%
        DataFrame
}







# Cellular barcode distributions ====
# Histogram function updated to use stored barcode counts rather than by
# accessing the files directly on the cluster. This approach allows us to
# generate multiple types of plots more easily, such as the new violin plot
# method.
.plot_cb_histogram <- function(run) {
    barcodes <- run$barcodes
    message("Generating barcode histograms...")
    pbmclapply(seq_along(barcodes), function(a) {
        name <- names(run$barcodes[a])
        sample_name <- run$metadata[name, "sample_name"]
        title <- paste(sample_name, name, sep = " : ")
        bcs <- run$barcodes[a] %>%
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

.plot_cb_violin <- function(run) {
    title <- deparse(substitute(run))
    barcodes <- unlist_barcodes(run)
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
#' @param run [bcbioSCDataSet].

#' @export
plot_barcodes <- function(run) {
    .plot_cb_violin(run) %>% show
    .plot_cb_histogram(run) %>% show
}



#' Unlist cellular barcodes
#'
#' Convert named list of cellular barcodes per sample to a data frame.
#'
#' @rdname unlist_barcodes
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param run [bcbioSCDataSet].
.unlist_barcodes <- function(run) {
    barcodes <- run$barcodes
    message("Converting nested barcodes to data frame...")
    pbmclapply(seq_along(barcodes), function(a) {
        barcodes[a] %>%
            as.data.frame %>%
            rownames_to_column %>%
            set_names(c("cellular_barcode", "reads")) %>%
            arrange(!!!syms(c("reads", "cellular_barcode"))) %>%
            mutate(log10_reads = log10(.data$reads),
                   sample_barcode = names(barcodes[a]))
    }) %>% bind_rows %>%
        left_join(run$metadata[, c("sample_barcode", "sample_name")],
                  by = "sample_barcode")
}
