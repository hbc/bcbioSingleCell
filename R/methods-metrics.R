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
        assert_are_identical(rownames(object), rownames(rowData))
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



.metrics.SE <- function(object) {  # nolint
    interestingGroups <- interestingGroups(object)
    sampleData <- sampleData(object)
    colData <- colData(object)
    colData[["cellID"]] <- rownames(colData)
    sampleID <- DataFrame("sampleID" = cell2sample(object))
    colData <- cbind(colData, sampleID)
    data <- merge(
        x = colData,
        y = sampleData,
        by = "sampleID",
        all.x = TRUE
    )
    data %>%
        as.data.frame() %>%
        .tidyMetrics() %>%
        # Ensure the metrics (colData) columns appear first
        .[, unique(c(colnames(colData), colnames(.)))] %>%
        column_to_rownames("cellID")
}



.tidyMetrics <- function(object) {
    assert_is_data.frame(object)
    # Ensure that rownames are set as a column before performing this chain
    object %>%
        # Enforce count columns as integers (e.g. `nUMI`)
        mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer) %>%
        # Coerce character vectors to factors, and drop levels
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels)
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
    signature("SingleCellExperiment"),
    .metrics.SE
)



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    function(object) {
        data <- .metrics.SE(object)
        # Ensure ident column is added
        assert_are_disjoint_sets(colnames(data), "ident")
        ident <- slot(object, "ident")
        cbind(data, ident)
    }
)
