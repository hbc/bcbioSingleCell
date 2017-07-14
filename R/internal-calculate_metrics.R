#' Calculate cellular barcode metrics summary
#'
#' This should only be performed during the initial run loading.
#'
#' @rdname calculate_metrics
#' @keywords internal
#'
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @param sparse_counts Sparse counts matrix.
#' @param annotable Annotable.
#' @param prefilter Whether to apply pre-filtering to the cellular barcodes.
#'
#' @return [matrix].
.calculate_metrics <- function(sparse_counts, annotable, prefilter = TRUE) {
    message("Calculating barcode metrics")
    cb_input <- ncol(sparse_counts)
    message(paste(cb_input, "cellular barcodes detected"))

    # Check for [Matrix::colSums()] methods support
    if (!"colSums,dgCMatrix-method" %in% methods(colSums)) {
        stop("dgCMatrix not supported in `colSums()`")
    }

    # Pull coding and mitochondrial genes from annotable
    coding_genes <- annotable %>%
        filter(.data[["broad_class"]] == "coding") %>%
        pull("ensgene")
    mito_genes <- annotable %>%
        filter(.data[["broad_class"]] == "mito") %>%
        pull("ensgene")

    metrics <- tibble(
        cellular_barcode = colnames(sparse_counts),
        umi_counts = Matrix::colSums(sparse_counts),
        genes_detected = Matrix::colSums(sparse_counts > 0L),
        coding_counts = Matrix::colSums(
            sparse_counts[rownames(sparse_counts) %in% coding_genes, ]),
        mito_counts = Matrix::colSums(
            sparse_counts[rownames(sparse_counts) %in% mito_genes, ])) %>%
        mutate(log10_genes_per_umi =
                   log10(.data[["genes_detected"]]) /
                   log10(.data[["umi_counts"]]),
               mito_ratio =
                   .data[["mito_counts"]] /
                   .data[["umi_counts"]])

    # Apply cellular barcode pre-filtering, if desired
    if (isTRUE(prefilter)) {
        metrics <- metrics %>%
            filter(.data[["umi_counts"]] >= 500L,
                   .data[["genes_detected"]] > 0L,
                   .data[["coding_counts"]] > 0L,
                   !is.na(.data[["log10_genes_per_umi"]]))
        message(paste(nrow(metrics), "cellular barcodes passed pre-filtering"))
    }

    metrics %>%
        as.data.frame %>%
        column_to_rownames %>%
        as.matrix
}
