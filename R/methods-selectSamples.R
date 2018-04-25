#' Select Samples
#'
#' @details Internally, pattern matching against sample and file names is
#'   applied using logical grep matching.
#'
#' @note Bracket based subsetting with `[` also works on `bcbioSingleCell`
#'   objects. In this case, provide cellular barcode identifiers for columns
#'   and Ensembl gene identifiers for rows.
#'
#' @name selectSamples
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase selectSamples
#'
#' @inheritParams general
#' @param ... Columns to use for grep pattern matching. Supply a named character
#'   vector containing the column name and the grep pattern.
#'
#' @return `bcbioSingleCell`.
#'
#' @seealso [sampleData()].
#'
#' @examples
#' # bcbioSingleCell ====
#' colnames(sampleData(bcb_small))
#' sampleName <- sampleData(bcb_small)[1L, "sampleName"]
#' print(sampleName)
#' selectSamples(bcb_small, sampleName = sampleName)
NULL



# Methods ======================================================================
#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("bcbioSingleCell"),
    function(object, ...) {
        metadata(object)[["selectSamples"]] <- TRUE

        # Here the `args` are captured as a named character vector. The
        # names of the arguments represent the column names. The value of the
        # arguments should be a string that can be used for logical grep
        # matching here internally.
        args <- list(...)
        checkAtomic <- vapply(
            X = args,
            FUN = is.atomic,
            FUN.VALUE = logical(1L)
        )
        if (!all(isTRUE(as.logical(checkAtomic)))) {
            stop("Arguments must be atomic vectors")
        }

        # Match the arguments against the sample metadata
        sampleData <- sampleData(object)
        sampleData[["sampleID"]] <- rownames(sampleData)

        list <- mapply(
            col = names(args),
            value = args,
            function(col, value) {
            # Check that column is present
            if (!col %in% colnames(sampleData)) {
                stop(paste(col, "isn't present in metadata colnames"))
            }
            # Check that all items in argument are present
            if (!all(value %in% sampleData[[col]])) {
                missing <- value[which(!value %in% sampleData[[col]])]
                stop(paste(
                    deparse(col),
                    "metadata column doesn't contain",
                    toString(missing)
                ))
            }
            sampleData %>%
                .[.[[col]] == value, , drop = FALSE] %>%
                rownames()
        })
        sampleIDs <- Reduce(f = intersect, x = list)

        # Output to the user which samples matched, using the `sampleName`
        # metadata column, which is more descriptive than `sampleID`
        sampleNames <- sampleData %>%
            .[sampleIDs, "sampleName", drop = TRUE] %>%
            as.character() %>%
            sort() %>%
            unique()

        message(paste(
            length(sampleNames), "sample(s) matched:",
            toString(sampleNames)
        ))

        # Use the metrics data.frame to match the cellular barcodes
        metrics <- metrics(object) %>%
            .[.[["sampleID"]] %in% sampleIDs, , drop = FALSE]
        if (!nrow(metrics)) {
            stop("Failed to match metrics")
        }
        message(paste(nrow(metrics), "cells matched"))
        cells <- rownames(metrics)

        object[, cells]
    }
)
