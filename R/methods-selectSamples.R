#' Select Samples
#'
#' @rdname selectSamples
#' @name selectSamples
#'
#' @details Internally, pattern matching against sample and file names is
#'   applied using [str_detect()].
#'
#' @note Bracket based subsetting with `[` also works on [bcbioSingleCell]
#'   objects. In this case, provide cellular barcode identifiers for columns
#'   and Ensembl gene identifiers for rows.
#'
#' @inheritParams AllGenerics
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
.selectSamples <- function(object, ...) {
    sampleMetadata <- sampleMetadata(object)
    patterns <- list(...)
    list <- lapply(seq_along(patterns), function(a) {
        col <- names(patterns)[[a]]
        pattern <- patterns[[a]]
        if (is_string(pattern)) {
            # Use grep pattern matching on string
            match <- sampleMetadata %>%
                .[str_detect(.[[col]], pattern), , drop = FALSE]
        } else if (is.character(pattern)) {
            # Use exact matching if vector supplied
            match <- sampleMetadata %>%
                .[.[[col]] %in% pattern, , drop = FALSE]
        } else {
            stop("Selection argument must be a character")
        }
        match %>%
            pull("sampleID") %>%
            unique() %>%
            sort()
    })
    # Use `base::Reduce()` explicitly here instead? Warning about init missing
    sampleIDs <- Reduce(f = intersect, x = list) %>%
        sort()
    if (!length(sampleIDs)) {
        stop("No samples matched")
    }
    message(paste(length(sampleIDs), "samples matched:", toString(sampleIDs)))
    sampleMetadata <- sampleMetadata %>%
        .[.[["sampleID"]] %in% sampleIDs, , drop = FALSE]

    # Match the sample ID prefix in the cellular barcode columns of the matrix.
    # Here `cb` is short for "cellular barcodes".
    cellularBarcodePattern <- sampleIDs %>%
        paste0(collapse = "|") %>%
        paste0("^(", ., ")_")
    cellularBarcodeMatches <- colnames(object) %>%
        .[str_detect(., cellularBarcodePattern)]
    if (!length(cellularBarcodeMatches)) stop("No cellular barcodes matched")
    message(paste(length(cellularBarcodeMatches), "cellular barcodes"))

    # Return the bcbioSingleCell object
    sparseCounts <- assay(object) %>%
        .[, cellularBarcodeMatches]
    colData <- colData(object) %>%
        .[cellularBarcodeMatches, ]
    rowData <- rowData(object) %>%
        set_rownames(rownames(object))
    metadata <- metadata(object)
    metadata[["sampleMetadata"]] <- sampleMetadata
    # Stash that samples are a subset
    metadata[["subset"]] <- TRUE
    se <- prepareSummarizedExperiment(
        assays = list(sparseCounts),
        colData = colData,
        rowData = rowData,
        metadata = metadata)
    new("bcbioSingleCell", se)
    # This will drop the unfiltered cellular barcodes
}



# Methods ====
#' @rdname selectSamples
#' @export
setMethod("selectSamples", "bcbioSingleCellANY", .selectSamples)
