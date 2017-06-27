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
    cb_input <- ncol(sparse_counts)
    message(paste(cb_input, "cellular barcodes detected"))

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

    metrics <- data.frame(
        rowname = colnames(sparse_counts),
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
                   .data[["umi_counts"]]) %>%
        # Pre-filter low quality barcodes
        filter(.data[["umi_counts"]] >= 500L,
               .data[["genes_detected"]] > 0L,
               .data[["coding_counts"]] > 0L,
               !is.na(.data[["log10_genes_per_umi"]])) %>%
        column_to_rownames %>%
        as.matrix

    message(paste(nrow(metrics), "cellular barcodes passed pre-filtering"))
    metrics
}
