#' Raw reads per cellular barcode
#'
#' Read counts prior to UMI disambiguation.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @note Updated 2019-08-08.
#' @noRd
#'
#' @param list `list`.
#' Cellular barcodes per sample.
#'
#' @return `integer`.
#' Cell identifiers are the names and raw reads are the values.
.nRead <- function(list) {
    assert(
        is.list(list),
        hasNames(list),
        is.integer(list[[1L]]),
        hasNames(list[[1L]])
    )
    if (hasLength(list, n = 1L)) {
        list[[1L]]
    } else {
        ## This will unlist using a "." separator.
        ## Renaming "." to "_" in names.
        x <- unlist(list, use.names = TRUE)
        names(x) <- makeNames(names(x))
        x
    }
}



#' Obtain the raw, unfiltered cellular barcode read counts
#'
#' @note Updated 2021-09-03.
#' @noRd
#'
#' @return `DataFrame`.
.rawMetrics <- function(object) {
    assert(is(object, "bcbioSingleCell"))
    list <- metadata(object)[["cellularBarcodes"]]
    assert(
        is.list(list),
        msg = sprintf(
            fmt = paste0(
                "Object does not contain unfiltered cellular barcodes.\n",
                "Has '%s' been applied? This step drops them."
            ),
            "filterCells()"
        )
    )
    assert(
        is.list(list),
        hasNames(list)
    )
    list <- Map(
        sampleId = names(list),
        reads = list,
        f = function(sampleId, reads) {
            DataFrame(
                "sampleId" = as.factor(sampleId),
                "cellId" = as.factor(names(reads)),
                "nRead" = reads,
                row.names = NULL
            )
        }
    )
    data <- unlist(DataFrameList(list), use.names = FALSE)
    sampleData <- sampleData(object)
    sampleData[["sampleId"]] <- rownames(sampleData)
    data <- leftJoin(data, sampleData, by = "sampleId")
    assert(
        is(data, "DataFrame"),
        !hasRownames(data),
        isSubset(c("sampleId", "cellId", "nRead"), colnames(data)),
        is.integer(data[["nRead"]])
    )
    data
}
