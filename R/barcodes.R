#' Generate barcode metrics summary
#'
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param run [bcbioSCDataSet].
#'
#' @return Tibble grouped by sample name.
#' @export
barcode_metrics <- function(run) {
    # Check for [Matrix::colSums()]
    if (!"colSums,dgCMatrix-method" %in% methods(colSums)) {
        stop("Invalid colSums function loaded")
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

    tibble(
        # Unique: `file_name` + `sample_barcode` + `cellular_barcode`
        unique = colnames(counts),
        total_counts = colSums(counts),
        genes_detected = colSums(counts > 0),
        coding_counts = colSums(
            counts[rownames(counts) %in% coding, ]),
        mito_counts = colSums(
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
        tidy_select(.data$sample_name, everything()) %>%
        # Filter barcodes matching samples
        filter(!is.na(.data$sample_name)) %>%
        group_by(!!!syms(c("sample_name",  "sample_barcode"))) %>%
        arrange(desc(!!sym("total_counts")), .by_group = TRUE)
}



#' Filter cellular barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @author Michael Steinbaugh
#'
#' @param run [bcbioSCDataSet].
#' @param reads Minimum number of total read counts per cell.
#' @param min_genes Minimum number of genes detected.
#' @param max_genes Maximum number of genes detected.
#' @param mito_ratio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param novelty Minimum novelty score.
#' @param show Show summary statistics and plots.
#'
#' @return Filtered [bcbioSCDataSet].
#' @export
filter_cellular_barcodes <- function(
    run,
    reads = 1000,
    min_genes = 500,
    max_genes = 5000,
    mito_ratio = 0.2,
    novelty = 0.8,
    show = TRUE) {
    name <- deparse(substitute(run))
    run$metrics <- run$metrics %>%
        filter(.data$total_counts > !!reads,
               # [fix] include option to filter by `coding_counts`?
               .data$genes_detected > !!min_genes,
               .data$genes_detected < !!max_genes,
               .data$mito_ratio < !!mito_ratio,
               .data$log10_detected_per_count > !!novelty)
    run$filtered <- TRUE
    if (isTRUE(show)) {
        writeLines(c(
            paste(name, "filtering parameters:"),
            "",
            paste0("- `> ", reads, "` total read counts per cell"),
            paste0("- `> ", min_genes, "` genes per cell"),
            paste0("- `< ", max_genes, "` genes per cell"),
            paste0("- `< ", mito_ratio, "` mitochondrial abundance ratio"),
            paste0("- `> ", novelty, "` novelty score")))
        plot_cell_counts(run)
        plot_read_counts(run, min = reads)
        plot_genes_detected(run, min = min_genes, max = max_genes)
        plot_reads_vs_genes(run) %>% show
        plot_mito_ratio(run, max = mito_ratio)
        plot_novelty(run, min = novelty)
    }
    run
}



# Cellular barcode distributions ====
# Histogram function updated to use stored barcode counts rather than by
# accessing the files directly on the cluster. This approach allows us to
# generate multiple types of plots more easily, such as the new violin plot
# method.
plot_cb_histogram <- function(run) {
    barcodes <- run$barcodes
    pbmclapply(seq_along(barcodes), function(a) {
        metadata <- run$metadata[sample_barcodes,
                                 c("sample_barcode", "sample_name")]
        plot <- function(file_name, sample_name) {
            bcs <- barcodes[a] %>%
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
                ggtitle(sample_name) +
                scale_x_log10(
                    breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(~10^.x))) +
                xlab("number of reads assigned to a cell") +
                ylab("proportion of cells")
        }
        show(plot(barcodes[a]))
    }) %>% invisible
}

plot_cb_violin <- function(run) {
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
#'
#' @return [show()] [ggplot].
#' @export
plot_barcodes <- function(run) {
    show(plot_cb_histogram)
    show(plot_cb_violin)
}



#' Top barcodes
#'
#' @param run [bcbioSCDataSet].
#' @param n Number of barcodes to return per sample.
#'
#' @export
top_barcodes <- function(run, n = 2) {
    run$metrics %>%
        top_n(n, !!sym("total_counts"))
}



#' Unlist cellular barcodes
#'
#' Convert named list of cellular barcodes per sample to a data frame.
#'
#' @keywords internal
#' @author Michael Steinbaugh
#'
#' @param run [bcbioSCDataSet].
#'
#' @export
unlist_barcodes <- function(run) {
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
