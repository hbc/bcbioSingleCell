#' Filter Cellular Barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filter
#' @author Michael Steinbaugh
#'
#' @param min_umis Minimum number of UMI disambiguated counts per cell.
#' @param min_genes Minimum number of genes detected.
#' @param max_genes Maximum number of genes detected.
#' @param max_mito_ratio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param min_novelty Minimum novelty score.
#' @param show_report Show summary statistics report and plots.
#'
#' @return [bcbioSCSubset].
#' @export
setMethod("filter", "bcbioSCDataSet", function(
    object,
    min_umis = 1000L,
    min_genes = 500L,
    max_genes = NULL,
    max_mito_ratio = 0.1,
    min_novelty = 0.8,
    show_report = TRUE) {
    sparse_counts <- assay(object)

    # Cellular barcode count
    ncol(sparse_counts) %>%
        paste("cellular barcodes in dataset") %>%
        message

    # Subset metrics
    metrics <- metrics(object) %>% as("tibble")

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

    # Convert back to data.frame
    metrics <- metrics %>%
        as.data.frame %>%
        column_to_rownames

    # Filter the sparse counts matrix with metrics
    sparse_size_original <- object.size(sparse_counts)
    sparse_counts <- sparse_counts[, rownames(metrics)]
    sparse_size_filtered <- object.size(sparse_counts)
    message(paste("Sparse matrix shrunk from",
                  format(sparse_size_original, units = "auto"),
                  "to",
                  format(sparse_size_filtered, units = "auto")))

    # colData ====
    col_data <- colData(object) %>% .[rownames(metrics), ]
    rm(metrics)

    # rowData ====
    row_data <- rowData(object) %>%
        set_rownames(rownames(object))

    # Metadata ====
    metadata <- SimpleList(
        filtering_criteria = c(
            min_umis = min_umis,
            min_genes = min_genes,
            max_genes = max_genes,
            max_mito_ratio = max_mito_ratio,
            min_novelty = min_novelty),
        pipeline = metadata(object)[["pipeline"]],
        source_name = deparse(substitute(object)),
        sample_metadata = metadata(object)[["sample_metadata"]],
        interesting_groups = metadata(object)[["interesting_groups"]],
        ensembl_version = metadata(object)[["ensembl_version"]],
        genome_build = metadata(object)[["genome_build"]],
        annotable = metadata(object)[["annotable"]])

    # SummarizedExperiment ====
    se <- packageSE(
        assays = SimpleList(
            sparse_counts = sparse_counts),
        colData = col_data,
        rowData = row_data,
        metadata = metadata)
    object <- new("bcbioSCSubset", se)

    # Show summary statistics report and plots, if desired
    if (isTRUE(show_report)) {
        mdHeader("Filtering parameters", level = 2L)
        writeLines(c(
            paste("  - >=", min_umis, "UMI counts per cell"),
            paste("  - >=", min_genes, "genes per cell"),
            paste("  - <=", max_genes, "genes per cell"),
            paste("  - <=", max_mito_ratio, "mitochondrial abundance ratio"),
            paste("  - >=", min_novelty, "novelty score")))

        mdHeader("Filtered cell counts", level = 3L)
        show(plot_cell_counts(object))

        mdHeader("Filtered UMI counts per cell", level = 3L)
        show(plot_umis_per_cell(object))

        mdHeader("Filtered genes detected", level = 3L)
        show(plot_genes_per_cell(object))

        mdHeader("Filtered mitochondrial counts ratio", level = 3L)
        show(plot_mito_ratio(object))

        mdHeader("Filtered novelty score", level = 3L)
        show(plot_novelty(object))
    }

    object
})
