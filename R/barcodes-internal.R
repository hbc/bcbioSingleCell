#' Counts per cellular barcode
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param list `list`.
#'   Cellular barcodes per sample.
#'
#' @return `integer`.
#' Cell identifiers are the names and raw read counts are the values.
.nCount <- function(list) {
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



## Obtain the raw, unfiltered cellular barcode read counts ("nCount").
.rawMetrics <- function(object) {
    list <- metadata(object)[["cellularBarcodes"]]

    if (!hasLength(list)) {
        message("Object is filtered. Using `metrics()` return.")
        return(metrics(object))
    }

    assert(is.list(list), isNonEmpty(list))

    if (hasLength(list, n = 1L)) {
        data <- tibble(
            sampleID = names(list)[[1L]],
            cellID = names(list[[1L]]),
            nCount = list[[1L]]
        )
    } else {
        data <- mapply(
            x = list,
            sampleID = names(list),
            FUN = function(x, sampleID) {
                tibble(
                    sampleID = sampleID,
                    cellID = names(x),
                    nCount = x
                )
            },
            SIMPLIFY = FALSE,
            USE.NAMES = FALSE
        )
        data <- bind_rows(data)
    }

    sampleData <- sampleData(object) %>%
        as_tibble(rownames = "sampleID") %>%
        mutate_all(as.factor)

    data <- data %>%
        mutate(!!sym("sampleID") := as.factor(!!sym("sampleID"))) %>%
        left_join(sampleData, by = "sampleID") %>%
        group_by(!!sym("sampleID"))

    assert(
        isSubset(c("sampleID", "cellID", "nCount"), colnames(data)),
        is.integer(data[["nCount"]])

    )

    data
}
