#' Read Known Markers
#'
#' @rdname readKnownMarkers
#' @author Michael Steinbaugh
#'
#' @param object Gene markers file (CSV or Excel).
#' @param genomeBuild Genome build.
#' @param show Show [kable].
#'
#' @return [tibble].
#' @export
setMethod("readKnownMarkers", "character", function(
    object, genomeBuild, show = TRUE) {
    annotable <- annotable(genomeBuild)
    markers <- readFileByExtension(object) %>%
        camel %>%
        left_join(annotable, by = "symbol") %>%
        .[!is.na(.data[["ensgene"]]), ] %>%
        .[, c("cell_type", "symbol")] %>%
        distinct
    if (isTRUE(show)) {
        kable(markers, caption = "Known markers") %>% show
    }
    markers
})
