#' Filter cellular barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filter
#' @author Michael Steinbaugh
#'
#' @param object Primary object.
#' @param min_umis Minimum number of UMI disambiguated counts per cell.
#' @param min_genes Minimum number of genes detected.
#' @param max_genes Maximum number of genes detected.
#' @param max_mito_ratio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param min_novelty Minimum novelty score.
#' @param show_report Show summary statistics report and plots.
#'
#' @return [SCSubset].
#' @export
setMethod("filter", "bcbioSCDataSet", function(
    object,
    min_umis = 1000L,
    min_genes = 500L,
    max_genes = NULL,
    max_mito_ratio = 0.1,
    min_novelty = 0.8,
    show_report = TRUE) {
    sparse_counts <- counts(object)

    # Cellular barcode count
    ncol(sparse_counts) %>%
        paste("cellular barcodes in dataset") %>%
        message

    # Subset metrics
    metrics <- metrics(object) %>% rownames_to_column

    # Apply filtering criteria
    if (!is.null(min_umis)) {
        metrics <- filter(metrics,
                          .data[["umi_counts"]] >= !!min_umis)
    }
    if (!is.null(min_genes)) {
        metrics <- filter(metrics,
                          .data[["genes_detected"]] >= !!min_genes)
    }
    if (!is.null(max_genes)) {
        metrics <- filter(metrics,
                          .data[["genes_detected"]] <= !!max_genes)
    }
    if (!is.null(max_mito_ratio)) {
        metrics <- filter(metrics,
                          .data[["mito_ratio"]] <= !!max_mito_ratio)
    }
    if (!is.null(min_novelty)) {
        metrics <- filter(metrics,
                          .data[["log10_genes_per_umi"]] >= !!min_novelty)
    }
    if (!nrow(metrics)) {
        stop("No cellular barcodes passed filtering")
    }
    message(paste(nrow(metrics), "cellular barcodes passed filtering"))
    metrics <- column_to_rownames(metrics)

    # Filter the sparse counts matrix with metrics
    sparse_size_original <- object.size(sparse_counts) %>%
        format(units = "auto")
    sparse_counts <- sparse_counts[, rownames(metrics)]
    sparse_size_filtered <- object.size(sparse_counts) %>%
        format(units = "auto")
    message(paste("Sparse matrix size shrunk from",
                  sparse_size_original, "to", sparse_size_filtered))

    # colData ====
    col_data <- colData(object) %>% .[rownames(metrics), ]

    # rowData ====
    row_data <- .row_data(object)

    # Metadata ====
    metadata <- SimpleList(
        source = deparse(substitute(object)),
        sample_metadata = metadata(object)[["sample_metadata"]],
        interesting_groups = metadata(object)[["interesting_groups"]],
        filtering_criteria = c(
            umis = umis,
            genes = genes,
            mito_ratio = mito_ratio,
            novelty = novelty),
        date = Sys.Date(),
        wd = getwd(),
        hpc = detect_hpc(),
        session_info = sessionInfo())

    # SummarizedExperiment ====
    se <- .summarized_experiment(
        assays = SimpleList(
            sparse_counts = sparse_counts),
        col_data = col_data,
        row_data = row_data,
        metadata = metadata)
    object <- new("SCSubset", se)

    # Show summary statistics report and plots, if desired
    if (isTRUE(show_report)) {
        writeLines(c(
            "Filtering parameters:",
            paste0("- `>= ", umis, "` UMI counts per cell"),
            paste0("- `>= ", genes, "` genes per cell"),
            paste0("- `<= ", mito_ratio, "` mitochondrial abundance ratio"),
            paste0("- `>= ", novelty, "` novelty score")))
        plot_cell_counts(object)
        plot_umis_per_cell(object, min = min_umis)
        plot_genes_detected(object, min = min_genes, max = max_genes)
        plot_mito_ratio(object, max = max_mito_ratio)
        plot_novelty(object, min = min_novelty)
    }

    object
})
