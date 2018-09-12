#' Sample Barcode Metrics
#'
#' @name metrics
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump metrics
#' @export
#'
#' @inheritParams general
#' @param recalculate `boolean`. Force recalculation, using primary [counts()]
#'   matrix.
#'
#' @return `data.frame` with cellular barcodes as rows. If `recalculate = TRUE`,
#'   returns the original object class (e.g. `SingleCellExperiment`), with
#'   updated metrics in [colData()].
#'
#' @examples
#' x <- metrics(indrops_small)
#' glimpse(x)
#'
#' # Recalculate mode
#' x <- metrics(indrops_small, recalculate = TRUE)
#' print(x)
NULL



.metrics <- function(
    object,
    rowRanges = NULL,
    prefilter = FALSE
) {
    assert_has_rows(object)
    assert_is_a_bool(prefilter)

    message("Calculating cellular barcode metrics")
    message(paste(ncol(object), "cells detected"))

    codingGenes <- character()
    mitoGenes <- character()

    if (length(rowRanges)) {
        assert_is_all_of(rowRanges, "GRanges")
        rowData <- as.data.frame(rowRanges)
        assertHasRownames(rowData)
        rowData <- rowData[rownames(object), , drop = FALSE]
        rowData <- as_tibble(rowData, rownames = "rowname")
        if ("broadClass" %in% colnames(rowData)) {
            # Drop rows with NA broad class
            rowData <- filter(rowData, !is.na(!!sym("broadClass")))
            # Coding genes
            codingGenes <- rowData %>%
                filter(!!sym("broadClass") == "coding") %>%
                pull("rowname")
            message(paste(length(codingGenes), "coding genes"))
            # Mitochondrial genes
            mitoGenes <- rowData %>%
                filter(!!sym("broadClass") == "mito") %>%
                pull("rowname")
            message(paste(length(mitoGenes), "mitochondrial genes"))
        }
    } else {
        warning(paste(
            "Calculating metrics without gene annotations.\n",
            "`rowRanges` is required to determine:",
            "`nCoding`, `nMito`, `mitoRatio`."
        ), call. = FALSE)
    }

    # Following the Seurat `seurat@meta.data` conventions
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

    # Return
    data %>%
        # Enforce count columns as integers (e.g. `nUMI`)
        mutate_if(grepl("^n[A-Z]", colnames(.)), as.integer) %>%
        # Coerce character vectors to factors, and drop levels
        mutate_if(is.character, as.factor) %>%
        mutate_if(is.factor, droplevels) %>%
        as.data.frame() %>%
        column_to_rownames()
}



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups = NULL,
        recalculate = FALSE
    ) {
        validObject(object)

        colData <- colData(object)

        # Calculate metrics, if necessary
        if (isTRUE(recalculate)) {
            # Bind metrics to colData
            metrics <- .metrics(
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
            return(object)
        }

        data <- as.data.frame(colData)
        if (!all(metricsCols %in% colnames(data))) {
            stop(paste(
                "Missing metrics:",
                setdiff(metricsCols, colnames(data)),
                "Run `metrics(object, recalculate = TRUE)`."
            ))
        }
        if (!"sampleName" %in% colnames(data)) {
            data[["sampleID"]] <- factor("unknown")
            data[["sampleName"]] <- factor("unknown")
        }
        interestingGroups <- matchInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )
        data <- uniteInterestingGroups(
            object = data,
            interestingGroups = interestingGroups
        )
        assert_is_subset(
            x = c("sampleName", "interestingGroups"),
            y = colnames(data)
        )

        data
    }
)
