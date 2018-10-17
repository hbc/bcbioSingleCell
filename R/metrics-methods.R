#' @rdname metrics
#' @importFrom basejump metrics
#' @export
metrics <- basejump::metrics



#' Metrics
#'
#' @name metrics
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit basejump::metrics
#'
#' @inheritParams general
#' @param recalculate `boolean`. Force recalculation, using primary [counts()]
#'   matrix.
#'
#' @examples
#' data(indrops_small)
#' x <- metrics(indrops_small, recalculate = TRUE)
#' print(x)
NULL



.metrics.matrix <- function(  # nolint
    object,
    rowRanges = NULL,
    prefilter = FALSE
) {
    assert_is_any_of(object, c("matrix", "sparseMatrix"))
    assert_has_rows(object)
    assert_is_any_of(rowRanges, c("GRanges", "NULL"))
    assert_is_a_bool(prefilter)

    message("Calculating cellular barcode metrics.")
    message(paste(ncol(object), "cells detected."))

    codingGenes <- character()
    mitoGenes <- character()

    rowRangesMessage <- paste(
        "`rowRanges` is required to calculate:",
        "nCoding, nMito, mitoRatio"
    )
    missingBiotype <- function() {
        message(paste(
            "Calculating metrics without biotype information.",
            rowRangesMessage,
            sep = "\n"
        ))
    }

    # Calculate nCoding and nMito, which requires annotations.
    if (length(rowRanges) > 0L) {
        assert_is_all_of(rowRanges, "GRanges")

        setdiff <- setdiff(rownames(object), names(rowRanges))
        if (has_length(setdiff)) {
            warning(paste(
                "Genes missing in rowRanges.",
                rowRangesMessage,
                printString(setdiff),
                sep = "\n"
            ))
            # Slot the rowRanges with empty ranges for these genes.
            # The same approach is used in `makeSummarizedExperiment()`.
            rowRanges <- suppressWarnings(c(
                rowRanges,
                emptyRanges(
                    names = setdiff,
                    mcolsNames = colnames(mcols(rowRanges))
                )
            ))
        }

        # Subset ranges to match matrix.
        assert_is_subset(rownames(object), names(rowRanges))
        rowRanges <- rowRanges[rownames(object)]
        rowData <- as(rowRanges, "tbl_df")
        if ("broadClass" %in% colnames(rowData)) {
            # Drop rows with NA broad class
            rowData <- filter(rowData, !is.na(!!sym("broadClass")))
            # Coding genes
            codingGenes <- rowData %>%
                filter(!!sym("broadClass") == "coding") %>%
                pull("rowname")
            message(paste(length(codingGenes), "coding genes."))
            # Mitochondrial genes
            mitoGenes <- rowData %>%
                filter(!!sym("broadClass") == "mito") %>%
                pull("rowname")
            message(paste(length(mitoGenes), "mitochondrial genes."))
        } else {
            missingBiotype()
        }
    } else {
        missingBiotype()
    }

    # Following the Seurat `seurat@meta.data` naming conventions.
    data <- tibble(
        rowname = colnames(object),
        nUMI = colSums(object),
        nGene = colSums(object > 0L),
        nCoding = colSums(object[codingGenes, , drop = FALSE]),
        nMito = colSums(object[mitoGenes, , drop = FALSE])
    ) %>%
        mutate(
            log10GenesPerUMI = log10(!!sym("nGene")) / log10(!!sym("nUMI")),
            mitoRatio = !!sym("nMito") / !!sym("nUMI")
        ) %>%
        # Ensure `n`-prefixed count columns are integer.
        mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer)

    # Apply low stringency cellular barcode pre-filtering.
    # This keeps only cellular barcodes with non-zero genes.
    if (isTRUE(prefilter)) {
        data <- data %>%
            filter(!is.na(UQ(sym("log10GenesPerUMI")))) %>%
            filter(!!sym("nUMI") > 0L) %>%
            filter(!!sym("nGene") > 0L)
        message(paste(
            nrow(data), "/", ncol(object),
            "cellular barcodes passed pre-filtering",
            paste0("(", percent(nrow(data) / ncol(object)), ")")
        ))
    }

    # Return as DataFrame, containing cell IDs in the rownames.
    data %>%
        # Enforce count columns as integers (e.g. `nUMI`).
        mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer) %>%
        # Coerce character vectors to factors, and drop levels.
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        as("DataFrame")
}



.metrics.SCE <-  # nolint
    function(object, recalculate = FALSE) {
        validObject(object)
        if (isTRUE(recalculate)) {
            colData <- colData(object)
            metrics <- .metrics.matrix(
                object = counts(object),
                rowRanges = rowRanges(object)
            )
            colData <- colData[
                ,
                setdiff(colnames(colData), colnames(metrics)),
                drop = FALSE
                ]
            colData <- cbind(metrics, colData)
            colData(object) <- colData
            object
        } else {
            metrics(as(object, "SingleCellExperiment"))
        }
    }



#' @rdname metrics
#' @export
setMethod(
    f = "metrics",
    signature = signature("bcbioSingleCell"),
    definition = .metrics.SCE
)



#' @rdname metrics
#' @export
setMethod(
    f = "metrics",
    signature = signature("CellRanger"),
    definition = .metrics.SCE
)
