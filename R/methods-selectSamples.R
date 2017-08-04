#' Select Samples
#'
#' @rdname selectSamples
#' @name selectSamples
#'
#' @details Internally, pattern matching against sample and file names is
#'   applied using [str_detect()].
#'
#' @param ... Columns to use for grep pattern matching. Supply a named character
#'   vector containing the column name and the grep pattern.
#'
#' @return [bcbioSCSubset].
NULL



# Constructors ====
.selectSamples <- function(object, ...) {
    sampleMetadata <- sampleMetadata(object)
    patterns <- c(...)
    list <- lapply(seq_along(patterns), function(a) {
        sampleMetadata %>%
            .[str_detect(.[[names(patterns)[[a]]]], patterns[[a]]), ] %>%
            .[, "sampleID"] %>%
            unique %>%
            sort
    })
    sampleIDs <- Reduce(intersect, list) %>% sort
    if (!length(sampleIDs)) stop("No samples matched")
    message(paste(length(sampleIDs), "samples matched:", toString(sampleIDs)))
    sampleMetadata <- sampleMetadata %>%
        .[.[["sampleID"]] %in% sampleIDs, ]

    # Match the sample ID prefix in the cellular barcode columns of the matrix.
    # Here `cb` is short for "cellular barcodes".
    cbPattern <- sampleIDs %>%
        paste0(collapse = "|") %>%
        paste0("^(", ., ")_")
    cbMatches <- colnames(object) %>%
        .[str_detect(., cbPattern)]
    if (!length(cbMatches)) stop("No cellular barcodes matched")
    message(paste(length(cbMatches), "cellular barcodes"))

    # Return the bcbioSCSubset object
    sparseCounts <- assay(object) %>%
        .[, cbMatches]
    colData <- colData(object) %>%
        .[cbMatches, ]
    rowData <- rowData(object) %>%
        set_rownames(rownames(object))
    metadata <- metadata(object)
    metadata[["sampleMetadata"]] <- sampleMetadata
    se <- packageSE(
        sparseCounts,
        colData = colData,
        rowData = rowData,
        metadata = metadata)
    new("bcbioSCSubset", se)
}



# Methods ====
#' @rdname selectSamples
#' @export
setMethod("selectSamples", "bcbioSCDataSet", .selectSamples)



#' @rdname selectSamples
#' @export
setMethod("selectSamples", "bcbioSCSubset", .selectSamples)
