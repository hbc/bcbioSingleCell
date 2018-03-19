#' Select Samples
#'
#' @name selectSamples
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase selectSamples
#'
#' @details Internally, pattern matching against sample and file names is
#'   applied using logical grep matching.
#'
#' @note Bracket based subsetting with `[` also works on `bcbioSingleCell`
#'   objects. In this case, provide cellular barcode identifiers for columns
#'   and Ensembl gene identifiers for rows.
#'
#' @inheritParams general
#' @param ... Columns to use for grep pattern matching. Supply a named character
#'   vector containing the column name and the grep pattern.
#'
#' @return `bcbioSingleCell`.
#'
#' @examples
#' # bcbioSingleCell ====
#' sampleData(bcb_small) %>% glimpse()
#' selectSamples(bcb_small, sampleName = "M1_seq_rep_1")
NULL



# Constructors =================================================================
#' @importFrom bcbioBase prepareSummarizedExperiment
#' @importFrom dplyr mutate_if
#' @importFrom magrittr set_rownames
.selectSamples <- function(object, ...) {
    object <- .applyFilterCutoffs(object)
    metadata(object)[["selectSamples"]] <- TRUE

    # Here the `arguments` are captured as a named character vector. The names
    # of the arguments represent the column names. The value of the arguments
    # should be a string that can be used for logical grep matching here
    # internally.
    arguments <- list(...)
    checkCharacter <- vapply(
        X = arguments,
        FUN = is.character,
        FUN.VALUE = logical(1L)
    )
    if (!all(isTRUE(as.logical(checkCharacter)))) {
        abort("Arguments must be character vectors")
    }

    # Match the arguments against the sample metadata
    sampleData <- sampleData(object)
    list <- lapply(seq_along(arguments), function(a) {
        column <- names(arguments)[[a]]
        # Check that column is present
        if (!column %in% colnames(sampleData)) {
            abort(paste(column, "isn't present in metadata colnames"))
        }
        argument <- arguments[[a]]
        # Check that all items in argument are present
        if (!all(argument %in% sampleData[[column]])) {
            missing <- argument[which(!argument %in% sampleData[[column]])]
            abort(paste(
                column,
                "metadata column doesn't contain",
                toString(missing)
            ))
        }
        sampleData %>%
            .[.[[column]] %in% argument, "sampleID", drop = TRUE]
    })
    sampleIDs <- Reduce(f = intersect, x = list)

    # Output to the user which samples matched, using the `sampleName` metadata
    # column, which is more descriptive than `sampleID`
    sampleNames <- sampleData %>%
        .[.[["sampleID"]] %in% sampleIDs, "sampleName", drop = TRUE] %>%
        as.character() %>%
        sort() %>%
        unique()
    inform(paste(
        length(sampleNames), "sample(s) matched:",
        toString(sampleNames)
    ))

    # Use the metrics data.frame to match the cellular barcodes
    metrics <- metrics(object) %>%
        .[.[["sampleID"]] %in% sampleIDs, , drop = FALSE]
    if (!nrow(metrics)) {
        abort("Failed to match metrics")
    }
    inform(paste(nrow(metrics), "cells matched"))
    cells <- rownames(metrics)

    object[, cells]
}



# Methods ======================================================================
#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("bcbioSingleCell"),
    .selectSamples
)
