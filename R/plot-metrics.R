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
    if (!is.null(metadata(bcb)[["lanes"]])) {
        plot <- metrics(bcb) %>%
            left_join(bcb$metadata,
                      by = c("sample_name", "sample_barcode")) %>%
            group_by(!!!syms(c("sample_name", "lane"))) %>%
            summarise(cells = n()) %>%
            ggplot(
                aes_(x = ~lane,
                     y = ~cells,
                     fill = ~sample_name)) +
            facet_wrap(~sample_name)
    } else {
        plot <- bcb$metrics %>%
            group_by(!!sym("sample_name")) %>%
            summarise(cells = n()) %>%
            ggplot(
                aes_(x = ~sample_name,
                     y = ~cells,
                     fill = ~sample_name))
    }
    plot +
        labs(title = "cell counts",
             x = "sample",
             y = "cell count") +
        geom_bar(stat = "identity") +
        geom_text(vjust = -0.5, aes_(label = ~cells)) +
        theme(
            axis.text.x = element_text(angle = 90L, hjust = 1L),
            legend.position = "none")
}

#' @rdname qc_plots_metrics
#' @export
plot_cell_counts <- function(bcb) {
    show(.plot_cell_counts_barplot(bcb))
}



# Read counts ====
# [TODO] Take out "total" in plot title
.plot_read_counts_boxplot <- function(bcb, min, type = "total") {
    if (!type %in% c("coding", "total")) {
        stop("Invalid counts column prefix")
    }
    name <- paste(type, "read counts")
    metrics <- bcb$metrics %>%
        rename(counts = !!sym(paste(type, "counts", sep = "_")))
    metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~counts,
                 fill = ~sample_name)) +
        labs(title = paste(name, "boxplot"),
             x = "sample",
             y = name) +
        geom_boxplot() +
        geom_hline(color = warn_color, yintercept = min) +
        geom_label(
            data = aggregate(counts ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(counts)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")) +
        scale_y_log10() +
        expand_limits(y = 1L) +
        theme(
            axis.text.x = element_text(angle = 90L, hjust = 1L),
            legend.position = "none")
}

.plot_read_counts_histogram <- function(bcb, min, type = "total") {
    if (!type %in% c("coding", "total")) {
        stop("Invalid counts column prefix")
    }
    name <- paste(type, "read counts")
    metrics <- bcb$metrics %>%
        rename(counts = !!sym(paste(type, "counts", sep = "_")))
    metrics %>%
        ggplot(
            aes_(x = ~counts,
                 fill = ~sample_name)) +
        labs(title = paste(name, "histogram"),
             x = name) +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color, xintercept = min) +
        expand_limits(x = 0) +
        scale_y_sqrt() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none")
}

#' @rdname qc_plots_metrics
#' @export
plot_read_counts <- function(bcb, min = 1000L) {
    show(.plot_read_counts_boxplot(bcb, min))
    show(.plot_read_counts_histogram(bcb, min))
    # [TODO] Add coding / total ratio plot?
}



# Genes detected ====
.plot_genes_detected_boxplot <- function(bcb, min, max) {
    plot <- bcb$metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~genes_detected,
                 fill = ~sample_name)) +
        labs(title = "genes detected boxplot",
             x = "sample",
             y = "genes per cell") +
        geom_boxplot() +
        geom_hline(color = warn_color, yintercept = min) +
        geom_label(
            data = aggregate(genes_detected ~ sample_name,
                             bcb$metrics,
                             median),
            aes_(label = ~round(genes_detected)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")) +
        expand_limits(y = 0) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none")

    # Show max genes cutoff, if set
    if (!is.null(max)) {
        plot <- plot + geom_hline(color = warn_color, yintercept = max)
    }

    plot
}

.plot_genes_detected_histogram <- function(bcb, min, max) {
    plot <- bcb$metrics %>%
        ggplot(
            aes_(x = ~genes_detected,
                 fill = ~sample_name)) +
        labs(title = "genes detected histogram",
             x = "genes per cell") +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color, xintercept = min) +
        expand_limits(x = 0) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none")

    # Show max genes cutoff, if set
    if (!is.null(max)) {
        plot <- plot + geom_vline(color = warn_color, xintercept = max)
    }

    plot
}

#' @rdname qc_plots_metrics
#' @export
plot_genes_detected <- function(bcb, min = 500, max = NULL) {
    show(.plot_genes_detected_boxplot(bcb, min, max))
    show(.plot_genes_detected_histogram(bcb, min, max))
}



# Read counts vs. detected genes ====
#' @rdname qc_plots_metrics
#' @param intgroup Interesting group used to define colors.
#' @export
plot_reads_vs_genes <- function(bcb, intgroup = NULL) {
    if (is.null(intgroup)) {
        intgroup <- bcb$intgroup[1]
    }
    bcb$metrics %>%
        ggplot(
            aes_(x = ~total_counts,
                 y = ~genes_detected,
                 # [fix] use [quo()] instead?
                 color = as.name(intgroup))) +
        labs(title = "total counts vs. genes detected",
             x = "counts per cell",
             y = "genes per cell") +
        facet_wrap(~sample_name) +
        geom_point() +
        geom_smooth(method = "loess") +
        scale_x_log10() +
        scale_y_log10() +
        # expand_limits(x = 1, y = 1) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none")
}



# Mitochondrial abundance ====
.plot_mito_ratio_boxplot <- function(bcb, max) {
    bcb$metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~mito_ratio,
                 fill = ~sample_name)
        ) +
        labs(title = "mitochondrial abundance boxplot",
             x = "sample",
             y = "relative mitochondrial abundance") +
        geom_boxplot() +
        geom_hline(color = warn_color, yintercept = max) +
        geom_label(
            data = aggregate(mito_ratio ~ sample_name,
                             bcb$metrics,
                             median),
            aes_(label = ~round(mito_ratio, digits = 2)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")
        ) +
        expand_limits(y = 0) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none")
}

.plot_mito_ratio_histogram <- function(bcb, max) {
    bcb$metrics %>%
        ggplot(
            aes_(x = ~mito_ratio,
                 fill = ~sample_name)) +
        labs(title = "mitochondrial abundance histogram",
             x = "relative mitochondrial abundance") +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color, xintercept = max) +
        xlim(0, 1) +
        scale_y_sqrt() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none")
}

.plot_mito_ratio_scatterplot <- function(bcb) {
    bcb$metrics %>%
        ggplot(
            aes_(x = ~coding_counts,
                 y = ~mito_counts,
                 color = ~sample_name)) +
        labs(title = "mitochondrial abundance scatterplot",
             x = "mitochondrial counts",
             y = "coding counts") +
        facet_wrap(~sample_name) +
        geom_point() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none")
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
    bcb$metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~log10_detected_per_count,
                 fill = ~sample_name)) +
        labs(title = "novelty boxplot",
             x = "sample",
             y = "log10 genes detected per count") +
        geom_boxplot() +
        geom_hline(color = warn_color, yintercept = min) +
        geom_label(
            data = aggregate(log10_detected_per_count ~ sample_name,
                             bcb$metrics,
                             median),
            aes_(label = ~round(log10_detected_per_count, digits = 2L)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")) +
        expand_limits(y = 0) +
        theme(
            axis.text.x = element_text(angle = 90L, hjust = 1L),
            legend.position = "none")
}

.plot_novelty_histogram <- function(bcb, min) {
    bcb$metrics %>%
        ggplot(
            aes_(x = ~log10_detected_per_count,
                 fill = ~sample_name)) +
        labs(title = "novelty histogram",
             x = "log10 genes detected per count") +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color, xintercept = min) +
        xlim(0L, 1L) +
        scale_y_sqrt() +
        theme(
            axis.text.x = element_text(angle = 90L, hjust = 1L),
            legend.position = "none")
}

#' @rdname qc_plots_metrics
#' @export
plot_novelty <- function(bcb, min = 0.8) {
    show(.plot_novelty_boxplot(bcb, min))
    show(.plot_novelty_histogram(bcb, min))
}
