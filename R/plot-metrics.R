#' Cell metrics quality control plots
#'
#' Novelty score means log10 genes detected per count.
#'
#' @rdname qc_plots_metrics
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param bcb [bcbioSCDataSet].
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#'
#' @return [ggplot].



# Cell counts ====
.plot_cell_counts_barplot <- function(bcb) {
    interesting_group <- interesting_groups(bcb)[[1L]]
    cell_counts <- metrics(bcb) %>%
        group_by(!!sym("sample_id")) %>%
        summarise(cells = n())
    sample_metadata <- sample_metadata(bcb)

    ggplot(
        cell_counts,
        aes_(x = ~sample_id,
             y = ~cells,
             fill = as.name(interesting_group))) +
        labs(title = "cell counts",
             x = "sample",
             y = "cell count") +
        geom_bar(stat = "identity") +
        geom_text(vjust = -0.5, aes_(label = ~cells)) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

#' @rdname qc_plots_metrics
#' @export
plot_cell_counts <- function(bcb) {
    show(.plot_cell_counts_barplot(bcb))
}



# Read counts ====
.plot_umis_per_cell_boxplot <- function(bcb, min) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]

    # Add interesting group to count aggregation data frame, for coloring
    meta <- sample_metadata(bcb) %>%
        .[, c("sample_id", interesting_group)]
    median_umis <- aggregate(umi_counts ~ sample_id, metrics, median) %>%
        left_join(meta, by = "sample_id") %>%
        mutate(umi_counts = round(.data[["umi_counts"]]))

    ggplot(
        metrics,
        aes_(x = ~sample_id,
             y = ~umi_counts,
             fill = as.name(interesting_group))) +
        labs(title = "umi counts boxplot",
             x = "sample",
             y = "umis per cell") +
        geom_boxplot() +
        geom_label(
            data = median_umis,
            aes_(label = ~umi_counts),
            alpha = 0.75,
            label.padding = unit(0.1, "lines"),
            show.legend = FALSE) +
        geom_hline(color = warn_color, yintercept = min) +
        scale_y_log10() +
        expand_limits(y = 1L) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

.plot_umis_per_cell_histogram <- function(bcb, min) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~umi_counts,
             fill = as.name(interesting_group))) +
        labs(title = "umi counts histogram",
             x = "umis per cell") +
        facet_wrap(~sample_id) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color, xintercept = min) +
        expand_limits(x = 0L) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

#' @rdname qc_plots_metrics
#' @export
plot_umis_per_cell <- function(bcb, min = 1000L) {
    show(.plot_umis_per_cell_boxplot(bcb, min))
    show(.plot_umis_per_cell_histogram(bcb, min))
}



# Genes detected ====
.plot_genes_detected_boxplot <- function(bcb, min, max) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]

    # Add interesting group to count aggregation data frame, for coloring
    meta <- sample_metadata(bcb) %>%
        .[, c("sample_id", interesting_group)]
    median_genes <- aggregate(genes_detected ~ sample_id, metrics, median) %>%
        left_join(meta, by = "sample_id") %>%
        mutate(genes_detected = round(.data[["genes_detected"]]))

    plot <- ggplot(
        metrics,
        aes_(x = ~sample_id,
             y = ~genes_detected,
             fill = as.name(interesting_group))) +
        labs(title = "genes detected boxplot",
             x = "sample",
             y = "genes per cell") +
        geom_boxplot() +
        geom_hline(color = warn_color, yintercept = min) +
        geom_label(
            data = median_genes,
            aes_(label = ~genes_detected),
            alpha = 0.75,
            label.padding = unit(0.1, "lines"),
            show.legend = FALSE) +
        expand_limits(y = 0L) +
        theme(
            axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Show max genes cutoff, if set
    if (!is.null(max)) {
        plot <- plot + geom_hline(color = warn_color, yintercept = max)
    }

    plot
}

.plot_genes_detected_histogram <- function(bcb, min, max) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]
    plot <- ggplot(
        metrics,
        aes_(x = ~genes_detected,
             fill = as.name(interesting_group))) +
        labs(title = "genes detected histogram",
             x = "genes per cell") +
        facet_wrap(~sample_id) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color, xintercept = min) +
        expand_limits(x = 0L) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))

    # Show max genes cutoff, if set
    if (!is.null(max)) {
        plot <- plot + geom_vline(color = warn_color, xintercept = max)
    }

    plot
}

#' @rdname qc_plots_metrics
#' @export
plot_genes_detected <- function(bcb, min = 500L, max = NULL) {
    show(.plot_genes_detected_boxplot(bcb, min, max))
    show(.plot_genes_detected_histogram(bcb, min, max))
}



# Read counts vs. detected genes ====
.plot_umis_vs_genes <- function(bcb) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~umi_counts,
             y = ~genes_detected,
             color = as.name(interesting_group))) +
        labs(title = "umis vs. genes detected",
             x = "umis per cell",
             y = "genes per cell") +
        facet_wrap(~sample_id) +
        geom_point() +
        geom_smooth(method = "loess") +
        scale_x_log10() +
        scale_y_log10() +
        # expand_limits(x = 1, y = 1) +
        theme(
            axis.text.x = element_text(angle = 90L, hjust = 1L),
            legend.position = "none")
}

#' @rdname qc_plots_metrics
#' @export
plot_umis_vs_genes <- function(bcb) {
    show(.plot_umis_vs_genes(bcb))
}



# Mitochondrial abundance ====
.plot_mito_ratio_boxplot <- function(bcb, max) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]

    # Add interesting group to count aggregation data frame, for coloring
    meta <- sample_metadata(bcb) %>%
        .[, c("sample_id", interesting_group)]
    median_mito_ratio <-
        aggregate(mito_ratio ~ sample_id, metrics, median) %>%
        left_join(meta, by = "sample_id") %>%
        mutate(mito_ratio = round(.data[["mito_ratio"]], digits = 3L))

    ggplot(
        metrics,
        aes_(x = ~sample_id,
             y = ~mito_ratio,
             fill = as.name(interesting_group))
    ) +
        labs(title = "mitochondrial abundance boxplot",
             x = "sample",
             y = "relative mitochondrial abundance") +
        geom_boxplot() +
        geom_hline(color = warn_color, yintercept = max) +
        geom_label(
            data = median_mito_ratio,
            aes_(label = ~mito_ratio),
            alpha = 0.75,
            label.padding = unit(0.1, "lines"),
            show.legend = FALSE) +
        expand_limits(y = 0L) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

.plot_mito_ratio_histogram <- function(bcb, max) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~mito_ratio,
             fill = as.name(interesting_group))) +
        labs(title = "mitochondrial abundance histogram",
             x = "relative mitochondrial abundance") +
        facet_wrap(~sample_id) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color, xintercept = max) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

.plot_mito_ratio_scatterplot <- function(bcb) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~coding_counts,
             y = ~mito_counts,
             color = as.name(interesting_group))) +
        labs(title = "mitochondrial abundance scatterplot",
             x = "mitochondrial counts",
             y = "coding counts") +
        facet_wrap(~sample_id) +
        geom_point() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

#' @rdname qc_plots_metrics
#' @export
plot_mito_ratio <- function(bcb, max = 0.2) {
    show(.plot_mito_ratio_boxplot(bcb, max))
    show(.plot_mito_ratio_histogram(bcb, max))
    show(.plot_mito_ratio_scatterplot(bcb))
}



# Novelty ====
.plot_novelty_boxplot <- function(bcb, min) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]

    # Add interesting group to count aggregation data frame, for coloring
    meta <- sample_metadata(bcb) %>%
        .[, c("sample_id", interesting_group)]
    median_novelty <-
        aggregate(log10_genes_per_umi ~ sample_id, metrics, median) %>%
        left_join(meta, by = "sample_id") %>%
        mutate(log10_genes_per_umi =
                   round(.data[["log10_genes_per_umi"]], digits = 3L))

    ggplot(
        metrics,
        aes_(x = ~sample_id,
             y = ~log10_genes_per_umi,
             fill = as.name(interesting_group))) +
        labs(title = "novelty boxplot",
             x = "sample",
             y = "log10 genes per umi (novelty score)") +
        geom_boxplot() +
        geom_hline(color = warn_color, yintercept = min) +
        geom_label(
            data = median_novelty,
            aes_(label = ~log10_genes_per_umi),
            alpha = 0.75,
            label.padding = unit(0.1, "lines"),
            show.legend = FALSE) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

.plot_novelty_histogram <- function(bcb, min) {
    metrics <- metrics(bcb)
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~log10_genes_per_umi,
             fill = as.name(interesting_group))) +
        labs(title = "novelty histogram",
             x = "log10 genes per umi (novelty score)") +
        facet_wrap(~sample_id) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color, xintercept = min) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

#' @rdname qc_plots_metrics
#' @export
plot_novelty <- function(bcb, min = 0.8) {
    show(.plot_novelty_boxplot(bcb, min))
    show(.plot_novelty_histogram(bcb, min))
}
