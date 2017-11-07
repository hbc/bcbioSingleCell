#' Calculate Cellular Barcode Metrics Summary
#'
#' @rdname calculateMetrics
#' @name calculateMetrics
#' @family QC Metrics Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param annotable Annotable.
#' @param prefilter Whether to apply pre-filtering to the cellular barcodes.
#'
#' @note This should only be performed during the initial run loading.
#'
#' @return [data.frame].
NULL



# Constructors ====
#' Calculate Metrics from Sparse Counts Matrix
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom basejump pct
#' @importFrom dplyr filter mutate pull
#' @importFrom Matrix colSums
#' @importFrom scales percent
#' @importFrom tibble column_to_rownames tibble
.calculateMetricsSparse <- function(
    object,
    annotable,
    prefilter = TRUE) {
    message("Calculating barcode metrics")
    message(paste(ncol(object), "cellular barcodes detected"))

    # Check that all genes are in annotable
    missing <- rownames(object) %>%
        .[!rownames(object) %in% annotable[["ensgene"]]]

    if (identical(length(missing), nrow(object))) {
        stop(paste(
            "No genes in the counts matrix matched the annotable.",
            "Check to ensure 'organism' is correct."
        ), call. = FALSE)
    }

    if (length(missing) > 0) {
        warning(paste(
            length(missing),
            "genes missing in annotable used to calculate metrics",
            paste0("(", percent(length(missing) / ncol(object)), ")")
        ), call. = FALSE)
    }

    # Obtain detected coding and mitochondrial genes, using annotable
    codingGenesDetected <- annotable %>%
        filter(.data[["broadClass"]] == "coding") %>%
        pull("ensgene") %>%
        .[. %in% rownames(object)]
    if (length(codingGenesDetected) == 0) {
        stop("No coding genes detected", call. = FALSE)
    }
    mitoGenesDetected <- annotable %>%
        filter(.data[["broadClass"]] == "mito") %>%
        pull("ensgene") %>%
        .[. %in% rownames(object)]
    if (length(mitoGenesDetected) == 0) {
        warning("No mitochondrial genes detected", call. = FALSE)
    } else {
        message(paste(
            length(mitoGenesDetected),
            "mitochondrial genes detected"
        ))
    }

    metrics <- tibble(
        rowname = colnames(object),
        # Follow the Seurat `seurat@data.info` conventions
        nUMI = Matrix::colSums(object),
        nGene = Matrix::colSums(object > 0),
        nCoding = Matrix::colSums(
            object[codingGenesDetected, , drop = FALSE]),
        nMito = Matrix::colSums(
            object[mitoGenesDetected, , drop = FALSE])) %>%
        mutate(log10GenesPerUMI =
                   log10(.data[["nGene"]]) / log10(.data[["nUMI"]]),
               # Using `nUMI` here like in Seurat example
               mitoRatio =
                   .data[["nMito"]] / .data[["nUMI"]])

    # Apply low stringency cellular barcode pre-filtering, if desired
    if (isTRUE(prefilter)) {
        # Here we're only keeping cellular barcodes that have a gene detected
        metrics <- metrics %>%
            .[.[["nUMI"]] > 0, ] %>%
            .[.[["nGene"]] > 0, ] %>%
            .[!is.na(.[["log10GenesPerUMI"]]), ]
        message(paste(
            nrow(metrics), "cellular barcodes passed pre-filtering",
            paste0("(", percent(nrow(metrics) / ncol(object)), ")")
        ))
    }

    metrics %>%
        as.data.frame() %>%
        column_to_rownames()
}



# Methods ====
#' @rdname calculateMetrics
#' @export
setMethod(
    "calculateMetrics",
    signature("dgCMatrix"),
    .calculateMetricsSparse)



#' @rdname calculateMetrics
#' @export
setMethod(
    "calculateMetrics",
    signature("bcbioSingleCell"),
    function(
        object,
        prefilter = TRUE) {
        message("Recalculating cellular barcode metrics")
        .calculateMetricsSparse(
            assay(object),
            annotable = metadata(object)[["annotable"]],
            prefilter = prefilter)
    })
