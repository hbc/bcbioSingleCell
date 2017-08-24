#' Sample Barcode Metrics
#'
#' @rdname metrics
#' @name metrics
#'
#' @return [data.frame] with cellular barcodes as rows.
NULL



# Constructors ====
.metrics <- function(object) {
    umiType <- metadata(object)[["umiType"]]
    meta <- sampleMetadata(object) %>%
        .[, unique(c(metaPriorityCols, interestingGroups(object)))]
    colData <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column
    if (umiType == "surecell") {
        match <- colData[["rowname"]] %>%
            str_match("^(.+)_([ACGT]{6}_[ACGT]{6}_[ACGT]{6})$")
    } else {
        match <- colData[["rowname"]] %>%
            str_match("^(.+)_([ACGT]{8}_[ACGT]{8})$")
    }
    match <- match %>%
        as.data.frame %>%
        set_colnames(c("rowname", "sampleID", "cellularBarcode"))
    left_join(colData, match, by = "rowname") %>%
        left_join(meta, by = "sampleID") %>%
        as.data.frame %>%
        column_to_rownames
}



# Methods ====
#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSCDataSet", .metrics)



#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSCFiltered", .metrics)



#' @rdname metrics
#' @export
setMethod("metrics", "seurat", function(object) {
    object@data.info %>% camel
})
