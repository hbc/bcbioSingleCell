#' Sample Barcode Metrics
#'
#' @name metrics
#' @family Data Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom bcbioBase metrics
#'
#' @inheritParams general
#' @param rowData Data describing the rows of the object.
#' @param prefilter Whether to apply pre-filtering to the cellular barcodes.
#'
#' @return `data.frame` with cellular barcodes as rows.
#'
#' @examples
#' # bcbioSingleCell ====
#' metrics(bcb_small) %>% glimpse()
#'
#' # seurat ====
#' metrics(seurat_small) %>% glimpse()
#' metrics(pbmc_small) %>% glimpse()
#'
#' # dgCMatrix ====
#' mat <- counts(bcb_small)
#' metrics(mat) %>% glimpse()
NULL



# Constructors =================================================================
.metrics.dgCMatrix <- function(  # nolint
    object,
    rowData = NULL,
    prefilter = TRUE
) {
    assert_has_rows(object)
    assert_is_a_bool(prefilter)

    inform("Calculating cellular barcode metrics")
    inform(paste(ncol(object), "cells detected"))

    if (!is.null(rowData)) {
        assertHasRownames(rowData)
        assert_are_intersecting_sets(
            x = rownames(object),
            y = rownames(rowData)
        )
        cols <- c("geneID", "broadClass")
        assert_is_subset(cols, colnames(rowData))
        rowData <- rowData %>%
            .[rownames(object), cols, drop = FALSE] %>%
            .[complete.cases(.), , drop = FALSE]
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
    } else {
        codingGenes <- character()
        mitoGenes <- character()
    }

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



.metrics.SCE <- function(object, interestingGroups) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    data <- colData(object)
    data <- uniteInterestingGroups(data, interestingGroups)
    as.data.frame(data)
}



# Methods ======================================================================
#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("dgCMatrix"),
    .metrics.dgCMatrix
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("bcbioSingleCell"),
    .metrics.SCE
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    .metrics.SCE
)
