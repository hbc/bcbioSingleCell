#' Cell metrics quality control plots.
#'
#' @rdname qc_plots_metrics
#'
#' @author Michael Steinbaugh
#' @author Rory Kirchner
#'
#' @param run bcbio-nextgen scRNA-seq run.
#'
#' @param min_genes Recommended minimum gene count cutoff.
#' @param max_genes Recommended maximum gene count cutoff.
#' @param mito_ratio Recommended maximum relative mitochondrial abundance cutoff
#'   value.
#' @param novelty Recommended novelty cutoff value.
#'
#' @param ... Passthrough parameters.
#'
#' @return ggplot2 object.



# Total cells ====

#' @rdname qc_plots_metrics
#' @description Total cells barplot.
#' @export
plot_total_cells <- function(run) {
    check_run(run)
    metrics <- run$metrics
    barplot <- metrics %>%
        group_by_(.dots = "sample_name") %>%
        summarize_(total_cells = ~n()) %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~total_cells,
                 fill = ~sample_name)
        ) +
        labs(title = "total number of cells",
             x = "sample name",
             y = "cell count (barcode cutoff applied)") +
        geom_bar(stat = "identity") +
        geom_text(vjust = -0.5,
                  aes_(label = ~total_cells)
        ) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(barplot)
}



# Total counts ====

#' @rdname qc_plots_metrics
#' @description Total counts histogram.
#' @export
plot_total_counts_histogram <- function(run) {
    check_run(run)
    metrics <- run$metrics
    histogram <- metrics %>%
        ggplot(
            aes_(x = ~total_counts,
                 fill = ~sample_name)
        ) +
        labs(title = "total RNA read counts histogram",
             x = "counts per cell") +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        expand_limits(x = 0) +
        scale_y_sqrt() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(histogram)
}

#' @rdname qc_plots_metrics
#' @description Total counts boxplot.
#' @export
plot_total_counts_boxplot <- function(run) {
    check_run(run)
    metrics <- run$metrics
    boxplot <- metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~total_counts,
                 fill = ~sample_name)
        ) +
        labs(title = "total RNA read counts boxplot",
             x = "sample name",
             y = "counts per cell") +
        geom_boxplot() +
        geom_label(
            data = aggregate(total_counts ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(total_counts)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")
        ) +
        scale_y_log10() +
        expand_limits(y = 1) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(boxplot)
}



# Genes detected ====

#' @rdname qc_plots_metrics
#' @description Genes detected boxplot.
#' @export
plot_genes_detected_boxplot <- function(
    run,
    min_genes = get("min_genes", envir = parent.frame()),
    max_genes = get("max_genes", envir = parent.frame())) {
    check_run(run)
    metrics <- run$metrics
    boxplot <- metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~genes_detected,
                 fill = ~sample_name)
        ) +
        labs(title = "genes detected boxplot",
             x = "sample name",
             y = "genes per cell") +
        geom_boxplot() +
        geom_hline(color = warn_color,
                   yintercept = min_genes) +
        geom_hline(color = warn_color,
                   yintercept = max_genes) +
        geom_label(
            data = aggregate(genes_detected ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(genes_detected)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")
        ) +
        expand_limits(y = 0) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(boxplot)
}

#' @rdname qc_plots_metrics
#' @description Genes detected histogram.
#' @export
plot_genes_detected_histogram <- function(
    run,
    min_genes = get("min_genes", envir = parent.frame()),
    max_genes = get("max_genes", envir = parent.frame())) {
    check_run(run)
    metrics <- run$metrics
    histogram <- metrics %>%
        ggplot(
            aes_(x = ~genes_detected,
                 fill = ~sample_name)
        ) +
        labs(title = "genes detected histogram",
             x = "genes per cell") +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color,
                   xintercept = min_genes) +
        geom_vline(color = warn_color,
                   xintercept = max_genes) +
        expand_limits(x = 0) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(histogram)
}



# Total vs. detected ====

#' @rdname qc_plots_metrics
#' @description Total counts vs. genes detected plot.
#'
#' @param colorby Column to color the points by.
#'
#' @export
plot_total_vs_detected <- function(run, colorby = "sample_name") {
    check_run(run)
    metrics <- run$metrics
    plot <- metrics %>%
        ggplot(
            aes_(x = ~total_counts,
                 y = ~genes_detected,
                 color = as.name(colorby))
        ) +
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
            legend.position = "none"
        )
    return(plot)
}



# Mitochondrial abundance ====

#' @rdname qc_plots_metrics
#' @description Mitochondrial abundance histogram.
#' @export
plot_mito_histogram <- function(
    run,
    mito_ratio = get("mito_ratio", envir = parent.frame())) {
    check_run(run)
    metrics <- run$metrics
    histogram <- metrics %>%
        ggplot(
            aes_(x = ~mito_ratio,
                 fill = ~sample_name)
        ) +
        labs(title = "mitochondrial gene abundance histogram",
             x = "relative mitochondrial abundance") +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color,
                   xintercept = mito_ratio) +
        xlim(0, 1) +
        scale_y_sqrt() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(histogram)
}

#' @rdname qc_plots_metrics
#' @description Mitochondrial abundance boxplot.
#' @export
plot_mito_boxplot <- function(
    run,
    mito_ratio = get("mito_ratio", envir = parent.frame())) {
    check_run(run)
    metrics <- run$metrics
    boxplot <- metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~mito_ratio,
                 fill = ~sample_name)
        ) +
        labs(title = "",
             x = "sample name",
             y = "relative mitochondrial abundance") +
        geom_boxplot() +
        geom_hline(color = warn_color,
                   yintercept = mito_ratio) +
        geom_label(
            data = aggregate(mito_ratio ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(mito_ratio, digits = 2)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")
        ) +
        expand_limits(y = 0) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(boxplot)
}

#' @rdname qc_plots_metrics
#' @description Mitochondrial abundance scatterplot.
#' @export
plot_mito_scatterplot <- function(run, mito_ratio = NULL) {
    check_run(run)
    metrics <- run$metrics
    scatterplot <- metrics %>%
        ggplot(
            aes_(x = ~coding_counts,
                 y = ~mito_counts,
                 color = ~sample_name)
        ) +
        labs(title = "mitochondrial gene abundance scatterplot",
             x = "counts in mitochondrial genes",
             y = "counts in coding genes") +
        facet_wrap(~sample_name) +
        geom_point() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(scatterplot)
}



# Novelty ====

#' @rdname qc_plots_metrics
#' @description Novelty histogram (log10 genes detected per count).
#' @export
plot_novelty_histogram <- function(
    run,
    novelty = get("novelty", envir = parent.frame())) {
    check_run(run)
    metrics <- run$metrics
    histogram <- metrics %>%
        ggplot(
            aes_(x = ~log10_detected_per_count,
                 fill = ~sample_name)
        ) +
        labs(title = "novelty histogram",
             x = "log10 genes detected per count") +
        facet_wrap(~sample_name) +
        geom_histogram(bins = bins) +
        geom_vline(color = warn_color, xintercept = novelty) +
        xlim(0, 1) +
        scale_y_sqrt() +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(histogram)
}

#' @rdname qc_plots_metrics
#' @description Novelty boxplot (log10 genes detected per count).
#' @export
plot_novelty_boxplot <- function(
    run,
    novelty = get("novelty", envir = parent.frame())) {
    check_run(run)
    metrics <- run$metrics
    boxplot <- metrics %>%
        ggplot(
            aes_(x = ~sample_name,
                 y = ~log10_detected_per_count,
                 fill = ~sample_name)
        ) +
        labs(title = "novelty boxplot",
             x = "sample name",
             y = "log10 genes detected per count") +
        geom_boxplot() +
        geom_hline(color = warn_color, yintercept = novelty) +
        geom_label(
            data = aggregate(log10_detected_per_count ~ sample_name,
                             metrics,
                             median),
            aes_(label = ~round(log10_detected_per_count, digits = 2)),
            alpha = 0.75,
            label.padding = unit(0.1, "lines")
        ) +
        expand_limits(y = 0) +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "none"
        )
    return(boxplot)
}



# RMarkdown/knit chunk wrappers ====

#' @rdname qc_plots_metrics
#' @description Plot total counts (RMarkdown chunk wrapper).
#' @export
plot_total_counts <- function(...) {
    show(plot_total_counts_histogram(...))
    show(plot_total_counts_boxplot(...))
}

#' @rdname qc_plots_metrics
#' @description Plot genes detected (RMarkdown chunk wrapper).
#' @export
plot_genes_detected <- function(...) {
    show(plot_genes_detected_histogram(...))
    show(plot_genes_detected_boxplot(...))
}

#' @rdname qc_plots_metrics
#' @description Plot mitochondrial counts (RMarkdown chunk wrapper).
#' @export
plot_mito <- function(...) {
    show(plot_mito_histogram(...))
    show(plot_mito_boxplot(...))
    show(plot_mito_scatterplot(...))
}

#' @rdname qc_plots_metrics
#' @description Plot novelty (RMarkdown chunk wrapper).
#' @export
plot_novelty <- function(...) {
    show(plot_novelty_histogram(...))
    show(plot_novelty_boxplot(...))
}
