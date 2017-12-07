#' Select Samples
#'
#' @rdname selectSamples
#' @name selectSamples
#'
#' @importFrom basejump selectSamples
#'
#' @details Internally, pattern matching against sample and file names is
#'   applied using logical grep matching.
#'
#' @note Bracket based subsetting with `[` also works on [bcbioSingleCell]
#'   objects. In this case, provide cellular barcode identifiers for columns
#'   and Ensembl gene identifiers for rows.
#'
#' @inheritParams AllGenerics
#'
#' @param ... Columns to use for grep pattern matching. Supply a named character
#'   vector containing the column name and the grep pattern.
#'
#' @return [bcbioSingleCell].
#'
#' @examples
#' \dontrun{
#' data(bcbFiltered)
#' # grep pattern matching with string
#' selectSamples(bcbFiltered, sampleName = "wt")
#'
#' # Exact name matching with character vector
#' selectSamples(bcbFiltered, sampleName = c("wt1", "wt2"))
#' }
NULL



# Constructors ====
#' @importFrom basejump prepareSummarizedExperiment
#' @importFrom dplyr mutate_if pull
#' @importFrom magrittr set_rownames
.selectSamples <- function(
    object,
    ...) {
    object <- .applyFilterCutoffs(object)
    metadata(object)[["selectSamples"]] <- TRUE

    # Here the `arguments` are captured as a named character vector. The names
    # of the arguments represent the column names. The value of the arguments
    # should be a string that can be used for logical grep matching here
    # internally.
    arguments <- list(...)
    checkCharacter <- vapply(arguments, is.character, FUN.VALUE = logical(1))
    if (!all(isTRUE(as.logical(checkCharacter)))) {
        stop("'Arguments must be character")
    }

    # Match the arguments against the sample metadata
    sampleMetadata <- sampleMetadata(object)
    list <- lapply(seq_along(arguments), function(a) {
        column <- names(arguments)[[a]]
        argument <- arguments[[a]]
        match <- sampleMetadata %>%
                .[.[[column]] %in% argument, , drop = FALSE]
        # Check for match failure
        if (!nrow(match)) {
            warning(paste(
                "Match failure:",
                paste(column, "=", argument)
            ), call. = FALSE)
            return(NULL)
        }
        pull(match, "sampleID")
    })
    sampleIDs <- Reduce(f = intersect, x = list)
    if (!length(sampleIDs)) {
        warning("No samples matched", call. = FALSE)
        return(NULL)
    }
    sampleIDs <- sort(unique(sampleIDs))

    # Filter the sample metadata data.frame to only contain matching samples
    sampleMetadata <- sampleMetadata %>%
        .[.[["sampleID"]] %in% sampleIDs, , drop = FALSE] %>%
        mutate_if(is.factor, droplevels)

    message(paste(
        length(sampleIDs), "sample(s) matched:",
        toString(sort(sampleMetadata[["sampleName"]]))
    ))

    # Use the metrics data.frame to match the cellular barcodes
    metrics <- metrics(object) %>%
        .[.[["sampleID"]] %in% sampleIDs, , drop = FALSE]
    if (!nrow(metrics)) {
        stop("Failed to match metrics", call. = FALSE)
    }
    message(paste(nrow(metrics), "cells matched"))
    cells <- rownames(metrics)

    object[, cells]
}



# Methods ====
#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("bcbioSingleCell"),
    .selectSamples)
