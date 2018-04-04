#' Calculate Cellular Barcode Metrics Summary
#'
#' @name calculateMetrics
#' @family Data Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#' @param rowData Data describing the rows of the object.
#' @param prefilter Whether to apply pre-filtering to the cellular barcodes.
#'
#' @note This step is performed automatically during the initial load.
#'
#' @return `data.frame`.
#'
#' @examples
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
    rowData <- as.data.frame(rowData)
    rownames(rowData) <- rowData[[1L]]
    # Allow resizing of the rowData here to support spike-ins
    assert_are_intersecting_sets(rownames(object), rownames(rowData))
    rowData <- rowData[rownames(object), , drop = FALSE]
    rownames(rowData) <- rownames(object)
    assert_is_a_bool(prefilter)

    inform("Calculating cellular barcode metrics")
    inform(paste(ncol(object), "cells detected"))

    codingGenes <- rowData %>%
        .[.[["broadClass"]] == "coding", "geneID", drop = TRUE] %>%
        na.omit()
    assert_all_are_non_missing_nor_empty_character(codingGenes)
    inform(paste(length(codingGenes), "coding genes detected"))

    mitoGenes <- rowData %>%
        .[.[["broadClass"]] == "mito", "geneID", drop = TRUE] %>%
        na.omit()
    assert_all_are_non_missing_nor_empty_character(mitoGenes)
    inform(paste(length(mitoGenes), "mitochondrial genes detected"))

    data <- tibble(
        "rowname" = colnames(object),
        # Follow the Seurat `seurat@data.info` conventions
        "nUMI" = Matrix::colSums(object),
        "nGene" = Matrix::colSums(object > 0L),
        "nCoding" = Matrix::colSums(
            object[codingGenes, , drop = FALSE]
        ),
        "nMito" = Matrix::colSums(
            object[mitoGenes, , drop = FALSE]
        )
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
        data <- data %>%
            .[.[["nUMI"]] > 0L, , drop = FALSE] %>%
            .[.[["nGene"]] > 0L, , drop = FALSE] %>%
            .[!is.na(.[["log10GenesPerUMI"]]), , drop = FALSE]
        inform(paste(
            nrow(data), "cellular barcodes passed pre-filtering",
            paste0("(", percent(nrow(data) / ncol(object)), ")")
        ))
    }

    data %>%
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
