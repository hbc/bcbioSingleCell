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
#'   Cellular barcodes per sample.
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
#' @note Updated 2019-08-12.
#' @noRd
#'
#' @return `tbl_df`.
.rawMetrics <- function(object) {
    assert(is(object, "bcbioSingleCell"))
    list <- metadata(object)[["cellularBarcodes"]]
    if (!is.list(list)) {
        stop(
            "Object does not contain unfiltered cellular barcodes.\n",
            "Has 'filterCells()' been applied? This step drops them."
        )
    }
    assert(
        is.list(list),
        isNonEmpty(list),
        hasNames(list)
    )
    list <- mapply(
        sampleID = names(list),
        reads = list,
        FUN = function(sampleID, reads) {
            DataFrame(
                sampleID = as.factor(sampleID),
                cellID = as.factor(names(reads)),
                nRead = reads,
                row.names = NULL
            )
        },
        SIMPLIFY = FALSE,
        USE.NAMES = TRUE
    )
    data <- unlist(DataFrameList(list), use.names = FALSE)
    sampleData <- sampleData(object)
    sampleData[["sampleID"]] <- rownames(sampleData)
    data <- left_join(data, sampleData, by = "sampleID")
    assert(
        is(data, "DataFrame"),
        !hasRownames(data),
        isSubset(c("sampleID", "cellID", "nRead"), colnames(data)),
        is.integer(data[["nRead"]])
    )
    as_tibble(data)
}
