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
#' @importFrom dplyr pull
#' @importFrom magrittr set_rownames
#' @importFrom rlang is_string
.selectSamples <- function(
    object,
    ...) {
    object <- .applyFilterCutoffs(object)

    # Here the `arguments` are captured as a named character vector. The names
    # of the arguments represent the column names. The value of the arguments
    # should be a string that can be used for logical grep matching here
    # internally.
    arguments <- list(...)
    sampleMetadata <- sampleMetadata(object)
    list <- lapply(seq_along(arguments), function(a) {
        column <- names(arguments)[[a]]
        pattern <- arguments[[a]]
        if (is_string(pattern)) {
            # Use grep pattern matching on string
            match <- sampleMetadata %>%
                .[grepl(x = .[[column]], pattern = pattern), , drop = FALSE]
        } else if (is.character(pattern)) {
            # Use exact matching if vector supplied
            match <- sampleMetadata %>%
                .[.[[column]] %in% pattern, , drop = FALSE]
        } else {
            stop("Selection argument must be a character", call. = FALSE)
        }
        match %>%
            pull("sampleID") %>%
            unique() %>%
            sort()
    })
    sampleIDs <- Reduce(f = intersect, x = list) %>%
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

    message(paste(length(cells), "cellular barcodes"))

    # Update the bcbio slot
    # Drop the unfiltered cellular barcodes. Only keep this in the main object
    # saved using `loadSingleCell()`.
    bcbio(bcb, "cellularBarcodes") <- NULL

    # Update the metadata slot
    metadata(bcb)[["allSamples"]] <- FALSE
    metadata(bcb)[["filterCells"]] <- cells
    metadata(bcb)[["filterGenes"]] <- genes
    metadata(bcb)[["sampleMetadata"]] <- sampleMetadata
    metadata(bcb)[["selectSamples"]] <- TRUE

    bcb
}



# Methods ====
#' @rdname selectSamples
#' @export
setMethod(
    "selectSamples",
    signature("bcbioSingleCell"),
    .selectSamples)
