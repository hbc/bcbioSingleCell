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
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' calculateMetrics(bcb) %>% glimpse()
#'
#' # dgCMatrix
#' counts <- counts(bcb)
#' annotable <- annotable(bcb)
#' calculateMetrics(counts, annotable = annotable) %>% glimpse()
NULL



# Constructors =================================================================
#' Calculate Metrics from Sparse Counts Matrix
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom dplyr filter mutate pull
#' @importFrom Matrix colSums
#' @importFrom scales percent
#' @importFrom tibble column_to_rownames tibble
.calculateMetricsSparse <- function(
    object,
    annotable = TRUE,
    prefilter = TRUE) {
    inform("Calculating barcode metrics")
    inform(paste(ncol(object), "cellular barcodes detected"))

    if (isTRUE(annotable)) {
        organism <- rownames(object) %>%
            .[[1]] %>%
            detectOrganism()
        annotable <- annotable(organism)
    }

    # Check that all genes are in annotable
    missing <- rownames(object) %>%
        .[!rownames(object) %in% annotable[["ensgene"]]]
    if (identical(length(missing), nrow(object))) {
        abort(paste(
            "No genes in the counts matrix matched the annotable.",
            "Check to ensure 'organism' is correct."
        ))
    }

    if (length(missing) > 0) {
        warn(paste(
            length(missing),
            "genes missing in annotable used to calculate metrics",
            paste0("(", percent(length(missing) / nrow(object)), ")")
        ))
    }

    # Obtain detected coding and mitochondrial genes, using annotable
    codingGenesDetected <- annotable %>%
        filter(.data[["broadClass"]] == "coding") %>%
        pull("ensgene") %>%
        .[. %in% rownames(object)]
    if (length(codingGenesDetected) == 0) {
        abort("No coding genes detected")
    }
    mitoGenesDetected <- annotable %>%
        filter(.data[["broadClass"]] == "mito") %>%
        pull("ensgene") %>%
        .[. %in% rownames(object)]
    if (length(mitoGenesDetected) == 0) {
        warn("No mitochondrial genes detected")
    } else {
        inform(paste(
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
        inform(paste(
            nrow(metrics), "cellular barcodes passed pre-filtering",
            paste0("(", percent(nrow(metrics) / ncol(object)), ")")
        ))
    }

    metrics %>%
        as.data.frame() %>%
        column_to_rownames()
}



# Methods ======================================================================
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
        inform("Recalculating cellular barcode metrics")
        .calculateMetricsSparse(
            assay(object),
            annotable = metadata(object)[["annotable"]],
            prefilter = prefilter)
    })
