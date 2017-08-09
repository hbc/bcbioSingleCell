#' Read Known Markers
#'
#' @rdname readMarkers
#' @name readMarkers
#'
#' @param object Gene markers file (CSV or Excel).
#' @param genomeBuild Genome build.
#' @param show Show [kable].
#'
#' @return [tibble].
NULL



# Methods ====
#' @rdname readMarkers
#' @export
setMethod("readMarkers", "character", function(
    object, genomeBuild, show = FALSE) {
    annotable <- annotable(genomeBuild)
    markers <- readFileByExtension(object) %>%
        camel
    # Ensure Ensemble gene identifiers are not set
    markers[["ensgene"]] <- NULL
    markers <- markers %>%
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
        .[, c("cellType", "symbol", "ensgene")] %>%
        arrange(!!!syms(c("cellType", "symbol"))) %>%
        distinct
    if (isTRUE(show)) {
        kable(markers, caption = "Known markers") %>% show
    }
    markers
})
