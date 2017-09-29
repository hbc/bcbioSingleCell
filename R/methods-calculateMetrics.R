#' Calculate Cellular Barcode Metrics Summary
#'
#' @rdname calculateMetrics
#' @name calculateMetrics
#' @family QC Metrics Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param annotable Annotable.
#' @param prefilter Whether to apply pre-filtering to the cellular barcodes.
#'
#' @note This should only be performed during the initial run loading.
#'
#' @return [matrix].
NULL



# Methods ====
#' @rdname calculateMetrics
#' @export
setMethod("calculateMetrics", "dgCMatrix", function(
    object, annotable, prefilter = TRUE) {
    message("Calculating barcode metrics")
    message(paste(ncol(object), "cellular barcodes detected"))

    # Check that all genes are in annotable
    missing <- rownames(object) %>%
        .[!rownames(object) %in% annotable[["ensgene"]]]
    if (length(missing) > 0) {
        warning(paste(
            length(missing), "genes missing in annotable",
            "used to calculate metrics:",
            toString(missing)
        ))
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
        # Follow the Seurat `seurat@data.info` conventions
        nUMI = Matrix::colSums(object),
        nGene = Matrix::colSums(object > 0),
        nCoding = Matrix::colSums(
            object[rownames(object) %in% codingGenes, ]),
        nMito = Matrix::colSums(
            object[rownames(object) %in% mitoGenes, ])) %>%
        mutate(log10GenesPerUMI =
                   log10(.data[["nGene"]]) / log10(.data[["nUMI"]]),
               mitoRatio =
                   .data[["nMito"]] / .data[["nUMI"]])

    # Apply low stringency cellular barcode pre-filtering, if desired
    if (isTRUE(prefilter)) {
        metrics <- metrics %>%
            .[.[["nUMI"]] > 0, ] %>%
            .[.[["nGene"]] > 0, ] %>%
            .[!is.na(.[["log10GenesPerUMI"]]), ]
        message(paste(nrow(metrics), "cellular barcodes passed pre-filtering"))
    }

    metrics %>%
        as.data.frame() %>%
        column_to_rownames() %>%
        as.matrix()
})
