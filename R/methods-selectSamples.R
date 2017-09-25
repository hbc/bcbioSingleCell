#' Select Samples
#'
#' @rdname selectSamples
#' @name selectSamples
#'
#' @details Internally, pattern matching against sample and file names is
#'   applied using [str_detect()].
#'
#' @inheritParams AllGenerics
#' @param ... Columns to use for grep pattern matching. Supply a named character
#'   vector containing the column name and the grep pattern.
#'
#' @return [bcbioSCFiltered].
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
    sampleIDs <- Reduce(intersect, list) %>% sort()
    if (!length(sampleIDs)) {
        stop("No samples matched")
    }
    message(paste(length(sampleIDs), "samples matched:", toString(sampleIDs)))
    sampleMetadata <- sampleMetadata %>%
        .[.[["sampleID"]] %in% sampleIDs, , drop = FALSE]

    # Match the sample ID prefix in the cellular barcode columns of the matrix.
    # Here `cb` is short for "cellular barcodes".
    cbPattern <- sampleIDs %>%
        paste0(collapse = "|") %>%
        paste0("^(", ., ")_")
    cbMatches <- colnames(object) %>%
        .[str_detect(., cbPattern)]
    if (!length(cbMatches)) stop("No cellular barcodes matched")
    message(paste(length(cbMatches), "cellular barcodes"))

    # Return the bcbioSCFiltered object
    sparseCounts <- assay(object) %>%
        .[, cbMatches]
    colData <- colData(object) %>%
        .[cbMatches, ]
    rowData <- rowData(object) %>%
        set_rownames(rownames(object))
    metadata <- metadata(object)
    metadata[["sampleMetadata"]] <- sampleMetadata
    se <- prepareSummarizedExperiment(
        sparseCounts,
        colData = colData,
        rowData = rowData,
        metadata = metadata)
    new("bcbioSCFiltered", se)
}



# Methods ====
#' @rdname selectSamples
#' @export
setMethod("selectSamples", "bcbioSCFiltered", .selectSamples)
