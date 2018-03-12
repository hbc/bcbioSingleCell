#' Calculate Cellular Barcode Metrics Summary
#'
#' @rdname calculateMetrics
#' @name calculateMetrics
#' @family QC Metrics Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param rowData Data describing the rows of the object.
#' @param prefilter Whether to apply pre-filtering to the cellular barcodes.
#'
#' @note This should only be performed during the initial run loading.
#'
#' @return [data.frame].
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' calculateMetrics(bcb) %>% glimpse()
#'
#' # dgCMatrix
#' counts <- counts(bcb)
#' rowData <- rowData(bcb)
#' calculateMetrics(counts, rowData = rowData) %>% glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr filter mutate mutate_if pull
#' @importFrom scales percent
#' @importFrom tibble column_to_rownames tibble
.calculateMetrics.dgCMatrix <- function(  # nolint
    object,
    rowData = TRUE,
    prefilter = TRUE) {
    inform("Calculating barcode metrics")
    inform(paste(ncol(object), "cellular barcodes detected"))

    if (isTRUE(rowData)) {
        organism <- rownames(object) %>%
            .[[1L]] %>%
            detectOrganism()
        rowData <- rowData(organism)
    }

    # Check that all genes are in rowData
    missing <- rownames(object) %>%
        .[!rownames(object) %in% rowData[["ensgene"]]]
    if (identical(length(missing), nrow(object))) {
        abort(paste(
            "No genes in the counts matrix matched the rowData.",
            "Check to ensure `organism` argument is correct."
        ))
    }

    if (length(missing) > 0L) {
        warn(paste(
            length(missing),
            "genes missing in rowData used to calculate metrics",
            paste0("(", percent(length(missing) / nrow(object)), ")")
        ))
    }

    # Obtain detected coding and mitochondrial genes, using rowData
    codingGenesDetected <- rowData %>%
        filter(.data[["broadClass"]] == "coding") %>%
        pull("ensgene") %>%
        .[. %in% rownames(object)]
    if (length(codingGenesDetected) == 0L) {
        abort("No coding genes detected")
    }
    mitoGenesDetected <- rowData %>%
        filter(.data[["broadClass"]] == "mito") %>%
        pull("ensgene") %>%
        .[. %in% rownames(object)]
    if (length(mitoGenesDetected) == 0L) {
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
        nGene = Matrix::colSums(object > 0L),
        nCoding = Matrix::colSums(
            object[codingGenesDetected, , drop = FALSE]
        ),
        nMito = Matrix::colSums(
            object[mitoGenesDetected, , drop = FALSE])
        ) %>%
        mutate(
            log10GenesPerUMI = log10(.data[["nGene"]]) /
                log10(.data[["nUMI"]]),
            # Using `nUMI` here like in Seurat example
            mitoRatio = .data[["nMito"]] / .data[["nUMI"]]
        ) %>%
        # Ensure count columns are integer. `colSums()` outputs as numeric.
        mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer)

    # Apply low stringency cellular barcode pre-filtering, if desired
    if (isTRUE(prefilter)) {
        # Here we're only keeping cellular barcodes that have a gene detected
        metrics <- metrics %>%
            .[.[["nUMI"]] > 0L, , drop = FALSE] %>%
            .[.[["nGene"]] > 0L, , drop = FALSE] %>%
            .[!is.na(.[["log10GenesPerUMI"]]), , drop = FALSE]
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
    .calculateMetrics.dgCMatrix
)



#' @rdname calculateMetrics
#' @export
setMethod(
    "calculateMetrics",
    signature("bcbioSingleCell"),
    function(
        object,
        prefilter = TRUE) {
        inform("Recalculating cellular barcode metrics")
        .calculateMetrics.dgCMatrix(
            assay(object),
            rowData = metadata(object)[["rowData"]],
            prefilter = prefilter
        )
    }
)
