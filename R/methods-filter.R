#' Filter cellular barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filter
#' @author Michael Steinbaugh
#'
#' @param object Primary object.
#' @param umis Minimum number of UMI disambiguated counts per cell.
#' @param genes Minimum number of genes detected.
#' @param mito_ratio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param novelty Minimum novelty score.
#' @param show_report Show summary statistics report and plots.
#'
#' @return [bcbioSCDataSet], with low quality cellular barcodes that
#'   don't pass quality control cutoffs removed.
#' @export
setMethod("filter", "bcbioSCDataSet", function(
    object,
    umis = 1000L,
    genes = 500L,
    mito_ratio = 0.1,
    novelty = 0.8,
    show_report = TRUE) {
    sparse_counts <- counts(object)

    # Cellular barcode count
    ncol(sparse_counts) %>%
        paste("cellular barcodes in dataset") %>%
        message

    # Subset metrics
    metrics <- metrics(object) %>%
        rownames_to_column %>%
        filter(.data[["umi_counts"]] >= !!umis,
               .data[["genes_detected"]] >= !!genes,
               .data[["mito_ratio"]] <= !!mito_ratio,
               .data[["log10_genes_per_umi"]] >= !!novelty) %>%
        column_to_rownames

    # Subset counts
    if (!nrow(metrics)) {
        stop("Barcode filtering was too stringent. No cells remaining.")
    }
    message(paste(nrow(metrics), "cellular barcodes passed filtering"))

    # Sparse
    sparse_size_original <- object.size(sparse_counts) %>%
        format(units = "auto")
    sparse_counts <- sparse_counts[, rownames(metrics)]
    sparse_size_filtered <- object.size(sparse_counts) %>%
        format(units = "auto")
    message(paste("Sparse matrix size shrunk from",
                  sparse_size_original, "to", sparse_size_filtered))

    # colData ====
    col_data <- colData(object) %>%
        .[rownames(metrics), ] %>%
        as("DataFrame")

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
    object <- .summarized_experiment(
        assays = SimpleList(
            sparse_counts = sparse_counts),
        col_data = col_data,
        row_data = row_data,
        metadata = metadata)

    # Show summary statistics report and plots, if desired
    if (isTRUE(show_report)) {
        writeLines(c(
            "Filtering parameters:",
            paste0("- `>= ", umis, "` UMI counts per cell"),
            paste0("- `>= ", genes, "` genes per cell"),
            paste0("- `<= ", mito_ratio, "` mitochondrial abundance ratio"),
            paste0("- `>= ", novelty, "` novelty score")))
        plot_cell_counts(object)
        plot_umis_per_cell(object, min = umis)
        plot_genes_detected(object, min = genes)
        plot_umis_vs_genes(object)
        plot_mito_ratio(object, max = mito_ratio)
        plot_novelty(object, min = novelty)
    }

    object
})
