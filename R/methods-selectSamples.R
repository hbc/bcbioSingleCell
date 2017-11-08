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
#' @importFrom rlang is_string
.selectSamples <- function(
    object,
    ...) {
    object <- .applyFilterCutoffs(object)
    cells <- colnames(object)
    genes <- rownames(object)

    # Here the `arguments` are captured as a named character vector. The names
    # of the arguments represent the column names. The value of the arguments
    # should be a string that can be used for logical grep matching here
    # internally.
    arguments <- list(...)
    checkCharacter <- vapply(arguments, is.character, FUN.VALUE = logical(1))
    checkCharacter <- as.logical(checkCharacter)
    if (!all(isTRUE(checkCharacter))) {
        stop("Arguments must be character", call. = FALSE)
    }

    # Convert all factors to strings for matching
    sampleMetadata <- sampleMetadata(object)

    list <- lapply(seq_along(arguments), function(a) {
        column <- names(arguments)[[a]]
        argument <- arguments[[a]]
        if (is_string(argument)) {
            # Use grep pattern matching on string
            match <- sampleMetadata %>%
                .[grepl(x = .[[column]],
                        pattern = argument,
                        ignore.case = FALSE),
                  , drop = FALSE]
        } else {
            # Use exact matching if vector supplied
            match <- sampleMetadata %>%
                .[.[[column]] %in% argument, , drop = FALSE]
        }
        # Check for match failure
        if (!nrow(match)) {
            stop(paste(
                "Match failure:",
                paste(column, "=", argument)
            ), call. = FALSE)
        }
        pull(match, "sampleID")
    })

    sampleIDs <- Reduce(f = intersect, x = list) %>%
        unique() %>%
        sort()
    if (!length(sampleIDs)) {
        stop("No samples matched", call. = FALSE)
    }

    # Filter the sample metadata data.frame to only contain matching samples
    sampleMetadata <- sampleMetadata %>%
        .[.[["sampleID"]] %in% sampleIDs, , drop = FALSE]

    message(paste(
        length(sampleIDs), "sample(s) matched:",
        toString(sort(sampleMetadata[["sampleName"]]))
    ))

    # Use the metrics data.frame to match the cellular barcodes
    metrics <- metrics(object, filterCells = TRUE) %>%
        .[.[["sampleID"]] %in% sampleIDs, , drop = FALSE]
    if (!nrow(metrics)) {
        stop("Failed to match metrics", call. = FALSE)
    }
    message(paste(nrow(metrics), "cellular barcodes"))
    cells <- rownames(metrics)

    # Now subset the object to contain only the cell matches
    subset <- object[, cells]

    # Update the bcbio slot
    cellularBarcodes <- bcbio(subset, "cellularBarcodes")
    if (!is.null(cellularBarcodes)) {
        bcbio(subset, "cellularBarcodes") <- cellularBarcodes[sampleID]
    }

    # Update the metadata slot
    metadata(subset)[["sampleMetadata"]] <- sampleMetadata
    metadata(subset)[["filterCells"]] <- cells
    metadata(subset)[["allSamples"]] <- FALSE
    metadata(subset)[["selectSamples"]] <- TRUE

    subset
}



# Methods ====
#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("bcbioSingleCell"),
    .selectSamples)
