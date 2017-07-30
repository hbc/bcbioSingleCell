#' Calculate Cellular Barcode Metrics Summary
#'
#' @rdname calculateMetrics
#'
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @param annotable Annotable.
#' @param prefilter Whether to apply pre-filtering to the cellular barcodes.
#'
#' @note This should only be performed during the initial run loading.
#'
#' @return [matrix].
#' @export
setMethod("calculateMetrics", "dgCMatrix", function(
    object, annotable, prefilter = TRUE) {
    message("Calculating barcode metrics")
    cbInput <- ncol(object)
    message(paste(cbInput, "cellular barcodes detected"))

    # Check that all genes are in annotable
    missingGenes <- rownames(object) %>%
        .[!rownames(object) %in% annotable[["ensgene"]]]
    if (length(missingGenes) > 0L) {
        warning(paste(length(missingGenes), "missing genes"))
    }

    # Check for [Matrix::colSums()] methods support
    if (!"colSums,dgCMatrix-method" %in% methods(colSums)) {
        stop("dgCMatrix not supported in `colSums()`", call. = FALSE)
    }

    # Extract coding and mitochondrial genes from annotable
    codingGenes <- annotable %>%
        .[.[["broadClass"]] == "coding", "ensgene"]
    mitoGenes <- annotable %>%
        .[.[["broadClass"]] == "mito", "ensgene"]

    metrics <- tibble(
        rowname = colnames(object),
        umiCounts = Matrix::colSums(object),
        genesDetected = Matrix::colSums(object > 0L),
        codingCounts = Matrix::colSums(
            object[rownames(object) %in% codingGenes, ]),
        mito_counts = Matrix::colSums(
            object[rownames(object) %in% mitoGenes, ])) %>%
        mutate(log10GenesPerUMI =
                   log10(.data[["genesDetected"]]) /
                   log10(.data[["umiCounts"]]),
               mitoRatio =
                   .data[["mitoCounts"]] /
                   .data[["umiCounts"]])

    # Apply low stringency cellular barcode pre-filtering, if desired
    if (isTRUE(prefilter)) {
        metrics <- metrics %>%
            .[.[["umiCounts"]] > 0L, ] %>%
            .[.[["genesDetected"]] > 0L, ] %>%
            .[.[["codingCounts"]] > 0L, ] %>%
            .[!is.na(.[["log10GenesPerUMI"]]), ]
        message(paste(nrow(metrics), "cellular barcodes passed pre-filtering"))
    }

    metrics %>%
        as.data.frame %>%
        column_to_rownames %>%
        as.matrix
})
