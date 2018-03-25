#' Calculate Cellular Barcode Metrics Summary
#'
#' @name calculateMetrics
#' @family Quality Control Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param rowData Data describing the rows of the object.
#' @param prefilter Whether to apply pre-filtering to the cellular barcodes.
#'
#' @note This should only be performed during the initial run loading.
#'
#' @return `data.frame`.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' calculateMetrics(bcb_small) %>% glimpse()
#'
#' # dgCMatrix ====
#' counts <- counts(bcb_small)
#' rowData <- rowData(bcb_small)
#' calculateMetrics(counts, rowData = rowData) %>% glimpse()
NULL



# Constructors =================================================================
.calculateMetrics.dgCMatrix <- function(  # nolint
    object,
    rowData,
    prefilter = TRUE
) {
    assert_is_any_of(rowData, c("data.frame", "DataFrame"))
    rowData <- as.data.frame(rowData)
    assert_are_intersecting_sets(rownames(object), rownames(rowData))
    assert_is_a_bool(prefilter)

    inform("Calculating cellular barcode metrics")
    inform(paste(ncol(object), "cells detected"))

    intersect <- intersect(rownames(object), rownames(rowData))
    setdiff <- setdiff(rownames(object), rownames(rowData))
    if (length(setdiff)) {
        warn(paste(
            length(setdiff),
            "genes missing in rowData",
            paste0("(", percent(length(setdiff) / nrow(object)), "):"),
            toString(setdiff)
        ))
    }
    rowData <- rowData[intersect, , drop = FALSE]

    # Obtain detected coding and mitochondrial genes, using rowData
    codingGenesDetected <- rowData %>%
        .[.[["broadClass"]] == "coding", , drop = FALSE] %>%
        pull("geneID")
    if (!length(codingGenesDetected)) {
        abort("No coding genes detected")
    }
    inform(paste(
        length(codingGenesDetected),
        "coding genes detected"
    ))
    mitoGenesDetected <- rowData %>%
        .[.[["broadClass"]] == "mito", , drop = FALSE] %>%
        pull("geneID")
    if (!length(mitoGenesDetected)) {
        abort("No mitochondrial genes detected")
    }
    inform(paste(
        length(mitoGenesDetected),
        "mitochondrial genes detected"
    ))

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
        prefilter = TRUE
    ) {
        inform("Recalculating cellular barcode metrics")
        calculateMetrics(
            object = assay(object),
            rowData = rowData(object),
            prefilter = prefilter
        )
    }
)
