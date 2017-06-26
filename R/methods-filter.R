# TODO Arrange plots into a grid?

#' Filter cellular barcodes
#'
#' Apply gene detection, mitochondrial abundance, and novelty score cutoffs to
#' cellular barcodes.
#'
#' @rdname filter
#' @author Michael Steinbaugh
#'
#' @param object [bcbioSCDataSet].
#' @param reads Minimum number of total read counts per cell.
#' @param min_genes Minimum number of genes detected.
#' @param max_genes Maximum number of genes detected.
#' @param mito_ratio Maximum relative mitochondrial abundance (`0-1` scale).
#' @param novelty Minimum novelty score.
#' @param show Show summary statistics and plots.
#'
#' @return [bcbioSCSubset].
#' @export
setMethod("filter", "bcbioSCDataSet", function(
    object,
    reads = 1000,
    min_genes = 500,
    max_genes = 5000,
    mito_ratio = 0.1,
    novelty = 0.8,
    show = TRUE) {
    name <- deparse(substitute(object))
    sparse_counts <- assay(object)

    # Cellular barcode count
    ncol(sparse_counts) %>%
        paste("cellular barcodes in dataset") %>%
        message

    # Subset metrics
    metrics <- metrics(object) %>%
        rownames_to_column %>%
        filter(.data[["total_counts"]] > !!reads,
               # FIXME Include option to filter by `coding_counts`?
               .data[["genes_detected"]] > !!min_genes,
               .data[["genes_detected"]] < !!max_genes,
               .data[["mito_ratio"]] < !!mito_ratio,
               .data[["log10_detected_per_count"]] > !!novelty) %>%
        column_to_rownames

    # Subset counts
    if (!nrow(metrics)) {
        stop("Barcode filtering was too stringent. No cells remaining.")
    }
    message(paste(nrow(metrics), "cellular barcodes passed filtering"))

    # ExpressionSet ====
    message("Packaging into ExpressionSet")
    # Sparse
    sparse_counts <- sparse_counts[, rownames(metrics)]
    sparse_size <- object.size(sparse_counts) %>% format(units = "auto")
    # Dense (currently required for assayData)
    dense_counts <- as.matrix(sparse_counts)
    dense_size <- object.size(dense_counts) %>% format(units = "auto")
    paste("Converting sparse to dense matrix:",
          sparse_size, "->", dense_size) %>%
        message
    pheno_data <- colData(object) %>%
        as.data.frame %>%
        .[rownames(metrics), ] %>%
        as("AnnotatedDataFrame")
    feature_data <- rowData(object) %>%
        as.data.frame %>%
        set_rownames(names(object)) %>%
        as("AnnotatedDataFrame")
    es <- ExpressionSet(
        assayData = dense_counts,
        phenoData = pheno_data,
        featureData = feature_data)
    bcb <- as(es, "bcbioSCSubset")

    # Unload counts from memory
    rm(sparse_counts, dense_counts)

    # Add additional slots
    filtering_criteria <- SimpleList(
        reads = reads,
        min_genes = min_genes,
        max_genes = max_genes,
        mito_ratio = mito_ratio,
        novelty = novelty)
    sample_metadata <- sample_metadata(object)
    interesting_groups <- metadata(object)[["interesting_groups"]]

    bcbio(bcb, "sample_metadata") <- sample_metadata
    bcbio(bcb, "interesting_groups") <- interesting_groups
    bcbio(bcb, "filtering_criteria") <- filtering_criteria

    message(paste(
        "bcbioSCSubset:",
        format(object.size(bcb), units = "auto")
    ))


    if (isTRUE(show)) {
        writeLines(c(
            paste(name, "filtering parameters:"),
            paste0("- `> ", reads, "` total read counts per cell"),
            paste0("- `> ", min_genes, "` genes per cell"),
            paste0("- `< ", max_genes, "` genes per cell"),
            paste0("- `< ", mito_ratio, "` mitochondrial abundance ratio"),
            paste0("- `> ", novelty, "` novelty score")))
        plot_cell_counts(bcb)
        plot_read_counts(bcb, min = reads)
        # FIXME non-numeric argument to binary operator?
        plot_genes_detected(bcb, min = min_genes, max = max_genes)
        plot_reads_vs_genes(bcb)
        plot_mito_ratio(bcb, max = mito_ratio)
        plot_novelty(bcb, min = novelty)
    }

    bcb
})
