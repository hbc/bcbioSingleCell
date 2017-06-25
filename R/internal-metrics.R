#' Generate barcode metrics summary
#'
#' @rdname metrics
#' @keywords internal
#' @author Rory Kirchner
#' @author Michael Steinbaugh
#'
#' @param sparse_counts Sparse counts matrix.
#' @param annotable Annotable.
#'
#' @return Tibble grouped by sample name.
.metrics <- function(sparse_counts, annotable) {
    message("Calculating barcode metrics")

    # Check for [Matrix::colSums()] methods support
    if (!"colSums,dgCMatrix-method" %in% methods(colSums)) {
        stop("dgCMatrix not supported in `colSums()`")
    }

    # Pull coding and mitochondrial genes from annotable
    annotable <- as.data.frame(annotable)
    coding_genes <- annotable %>%
        filter(.data[["broad_class"]] == "coding") %>%
        pull("ensgene")
    mito_genes <- annotable %>%
        filter(.data[["broad_class"]] == "mito") %>%
        pull("ensgene")

    data.frame(
        # Rowname: `file_name` + `sample_barcode` + `cellular_barcode`
        rowname = colnames(sparse_counts),
        total_counts = Matrix::colSums(sparse_counts),
        genes_detected = Matrix::colSums(sparse_counts > 0),
        coding_counts = Matrix::colSums(
            sparse_counts[rownames(sparse_counts) %in% coding_genes, ]),
        mito_counts = Matrix::colSums(
            sparse_counts[rownames(sparse_counts) %in% mito_genes, ])) %>%
        # Filter zero count barcodes
        filter(.data[["total_counts"]] > 0) %>%
        mutate(log10_detected_per_count =
                   log10(.data[["genes_detected"]]) /
                   log10(.data[["total_counts"]]),
               mito_ratio =
                   .data[["mito_counts"]] /
                   .data[["total_counts"]]) %>%
        column_to_rownames %>%
        as.matrix
}
