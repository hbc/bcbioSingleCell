#' Cell metrics quality control plots
#'
#' Novelty score means log10 genes detected per count.
#'
#' @rdname plot-metrics
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
    metrics <- metrics(bcb)
    cell_counts <- metrics %>%
        group_by(!!sym("sample_id")) %>%
        summarise(cells = n()) %>%
        left_join(sample_metadata(bcb), by = "sample_id")
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        cell_counts,
        aes_(x = ~sample_name,
             y = ~cells,
             fill = as.name(interesting_group))) +
        labs(title = "cell counts",
             x = "sample",
             y = "cell count") +
        facet_wrap(~file_name) +
        geom_bar(stat = "identity") +
        geom_text(vjust = -0.5, aes_(label = ~cells)) +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

#' @rdname plot-metrics
#' @export
plot_cell_counts <- function(bcb) {
    show(.plot_cell_counts_barplot(bcb))
}



# Read counts ====
.plot_umis_per_cell_boxplot <- function(bcb, min) {
    metrics <- metrics(bcb)
    median_umis <- aggregate(umi_counts ~ sample_id, metrics, median) %>%
        left_join(sample_metadata(bcb), by = "sample_id") %>%
        mutate(umi_counts = round(.data[["umi_counts"]]))
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~umi_counts,
             fill = as.name(interesting_group))) +
        labs(title = "umi counts boxplot",
             x = "sample",
             y = "umis per cell") +
        facet_wrap(~file_name) +
        geom_boxplot() +
        geom_label(data = median_umis,
                   aes_(label = ~umi_counts),
                   alpha = 0.75,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        geom_hline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   yintercept = min) +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

.plot_umis_per_cell_histogram <- function(bcb, min) {
    metrics <- metrics(bcb)
    ggplot(
        metrics,
        aes_(x = ~umi_counts,
             fill = ~sample_name)) +
        labs(title = "umi counts histogram",
             x = "umis per cell") +
        facet_wrap(~file_name) +
        geom_histogram(bins = bins) +
        geom_vline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   xintercept = min) +
        scale_x_log10() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

#' @rdname plot-metrics
#' @export
plot_umis_per_cell <- function(bcb, min = 1000L) {
    show(.plot_umis_per_cell_boxplot(bcb, min))
    show(.plot_umis_per_cell_histogram(bcb, min))
}



# Genes detected ====
.plot_genes_detected_boxplot <- function(bcb, min) {
    metrics <- metrics(bcb)
    median_genes <- aggregate(genes_detected ~ sample_id, metrics, median) %>%
        left_join(sample_metadata(bcb), by = "sample_id")
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~genes_detected,
             fill = as.name(interesting_group))) +
        labs(title = "genes detected boxplot",
             x = "sample",
             y = "genes per cell") +
        facet_wrap(~file_name) +
        geom_boxplot() +
        geom_hline(alpha = 0.5,
                   color = warn_color,
                   size = 2,
                   yintercept = min) +
        geom_label(data = median_genes,
                   aes_(label = ~round(genes_detected)),
                   alpha = 0.75,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

.plot_genes_detected_histogram <- function(bcb, min) {
    metrics <- metrics(bcb)
    ggplot(
        metrics,
        aes_(x = ~genes_detected,
             fill = ~sample_name)) +
        labs(title = "genes detected histogram",
             x = "genes per cell") +
        facet_wrap(~file_name) +
        geom_histogram(bins = bins) +
        geom_vline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   xintercept = min) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

#' @rdname plot-metrics
#' @export
plot_genes_detected <- function(bcb, min = 500L) {
    show(.plot_genes_detected_boxplot(bcb, min))
    show(.plot_genes_detected_histogram(bcb, min))
}



# Read counts vs. detected genes ====
.plot_umis_vs_genes <- function(bcb) {
    metrics <- metrics(bcb)
    ggplot(
        metrics,
        aes_(x = ~umi_counts,
             y = ~genes_detected,
             color = ~sample_name)) +
        labs(title = "umis vs. genes detected",
             x = "umis per cell",
             y = "genes per cell") +
        facet_wrap(~file_name) +
        geom_smooth(method = "lm", se = FALSE) +
        scale_x_log10() +
        scale_y_log10() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

#' @rdname plot-metrics
#' @export
plot_umis_vs_genes <- function(bcb) {
    show(.plot_umis_vs_genes(bcb))
}



# Mitochondrial abundance ====
.plot_mito_ratio_boxplot <- function(bcb, max) {
    metrics <- metrics(bcb)
    median_mito_ratio <-
        aggregate(mito_ratio ~ sample_id, metrics, median) %>%
        left_join(sample_metadata(bcb), by = "sample_id")
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~mito_ratio,
             fill = as.name(interesting_group))) +
        labs(title = "mitochondrial abundance boxplot",
             x = "sample",
             y = "relative mitochondrial abundance") +
        facet_wrap(~file_name) +
        geom_boxplot() +
        geom_hline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   yintercept = max) +
        geom_label(data = median_mito_ratio,
                   aes_(label = ~round(mito_ratio, digits = 2L)),
                   alpha = 0.75,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

.plot_mito_ratio_histogram <- function(bcb, max) {
    metrics <- metrics(bcb)
    ggplot(
        metrics,
        aes_(x = ~mito_ratio,
             fill = ~sample_name)) +
        labs(title = "mitochondrial abundance histogram",
             x = "relative mitochondrial abundance") +
        facet_wrap(~file_name) +
        geom_histogram(bins = bins) +
        geom_vline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   xintercept = max) +
        scale_x_sqrt() +
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

#' @rdname plot-metrics
#' @export
plot_mito_ratio <- function(bcb, max = 0.2) {
    show(.plot_mito_ratio_boxplot(bcb, max))
    show(.plot_mito_ratio_histogram(bcb, max))
    show(.plot_mito_ratio_scatterplot(bcb))
}



# Novelty ====
.plot_novelty_boxplot <- function(bcb, min) {
    metrics <- metrics(bcb)
    median_novelty <-
        aggregate(log10_genes_per_umi ~ sample_id, metrics, median) %>%
        left_join(sample_metadata(bcb), by = "sample_id")
    interesting_group <- interesting_groups(bcb)[[1L]]
    ggplot(
        metrics,
        aes_(x = ~sample_name,
             y = ~log10_genes_per_umi,
             fill = as.name(interesting_group))) +
        labs(title = "novelty boxplot",
             x = "sample",
             y = "log10 genes per umi (novelty score)") +
        facet_wrap(~file_name) +
        geom_boxplot() +
        geom_hline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   yintercept = min) +
        geom_label(data = median_novelty,
                   aes_(label = ~round(log10_genes_per_umi, digits = 2L)),
                   alpha = 0.75,
                   label.padding = unit(0.1, "lines"),
                   show.legend = FALSE) +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

.plot_novelty_histogram <- function(bcb, min) {
    metrics <- metrics(bcb)
    ggplot(
        metrics,
        aes_(x = ~log10_genes_per_umi,
             fill = ~sample_name)) +
        labs(title = "novelty histogram",
             x = "log10 genes per umi (novelty score)") +
        facet_wrap(~file_name) +
        geom_histogram(bins = bins) +
        geom_vline(alpha = 0.5,
                   color = warn_color,
                   size = 2L,
                   xintercept = min) +
        scale_x_sqrt() +
        scale_y_sqrt() +
        theme(axis.text.x = element_text(angle = 90L, hjust = 1L))
}

#' @rdname plot-metrics
#' @export
plot_novelty <- function(bcb, min = 0.8) {
    show(.plot_novelty_boxplot(bcb, min))
    show(.plot_novelty_histogram(bcb, min))
}
