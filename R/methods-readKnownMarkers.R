#' Read Known Markers
#'
#' @rdname readKnownMarkers
#' @name readKnownMarkers
#'
#' @param object Gene markers file (CSV or Excel).
#' @param genomeBuild Genome build.
#' @param show Show [kable].
#'
#' @return [tibble].
NULL



# Methods ====
#' @rdname readKnownMarkers
#' @export
setMethod("readKnownMarkers", "character", function(
    object, genomeBuild, show = TRUE) {
    annotable <- annotable(genomeBuild)
    markers <- readFileByExtension(object) %>%
        camel %>%
        # Remove rows that don't contain a symbol to map
        .[!is.na(.[["symbol"]]), ] %>%
        left_join(annotable, by = "symbol")
    if (any(is.na(markers[["ensgene"]]))) {
        missing <- markers %>%
            .[is.na(.[["ensgene"]]), ] %>%
            pull("symbol") %>%
            sort %>%
            unique
        stop(paste("Unmapped symbols:", toString(missing)))
    }
    markers <- markers %>%
        .[!is.na(.[["ensgene"]]), ] %>%
        .[, c("cellType", "symbol")] %>%
        distinct
    if (isTRUE(show)) {
        kable(markers, caption = "Known markers") %>% show
    }
    markers
})
