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
#' @importFrom basejump selectSamples
#' @export
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
#' object <- indrops_small
#' sample <- sampleNames(object) %>% head(1L)
#' print(sample)
#' selectSamples(object, sampleName = sample)
NULL



#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("SingleCellExperiment"),
    function(
        object,
        ...,
        prefilter = TRUE
    ) {
        validObject(object)
        assert_is_a_bool(prefilter)

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
        # Include sampleID for looping in other functions
        sampleData[["sampleID"]] <- rownames(sampleData)

        list <- mapply(
            col = names(args),
            arg = args,
            function(col, arg) {
                # Check that column is present
                if (!col %in% colnames(sampleData)) {
                    stop(paste(col, "isn't present in metadata colnames"))
                }
                # Check that all items in argument are present
                if (!all(arg %in% sampleData[[col]])) {
                    missing <- arg[which(!arg %in% sampleData[[col]])]
                    stop(paste(
                        deparse(col),
                        "metadata column doesn't contain",
                        toString(missing)
                    ))
                }
                sampleData %>%
                    .[.[[col]] %in% arg, , drop = FALSE] %>%
                    rownames()
            },
            SIMPLIFY = FALSE,
            USE.NAMES = TRUE
        )
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

        # Use the metrics `data.frame` to match the cellular barcodes
        metrics <- metrics(object)
        assert_is_subset("sampleID", colnames(metrics))
        cells <- metrics %>%
            rownames_to_column("cellID") %>%
            filter(!!sym("sampleID") %in% !!sampleIDs) %>%
            pull("cellID")
        message(paste(length(cells), "cells matched"))

        # Note that this step will drop the raw cellular barcodes list
        object <- object[, cells]

        # Run low-stringency filtering to drop zero-count genes
        if (isTRUE(prefilter)) {
            object <- filterCells(object)
        }

        object
    }
)
