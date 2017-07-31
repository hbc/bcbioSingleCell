#' Top Markers
#'
#' @rdname topMarkers
#'
#' @param n Number of genes per cluster.
#' @param show Show [kable].
#'
#' @return [tibble].
#' @export
setMethod("topMarkers", "data.frame", function(object, n = 4L, show = TRUE) {
    markers <- .groupMarkers(object) %>%
        top_n(n = n, wt = .data[["avgDiff"]])
    if (isTRUE(show)) {
        kable(markers,
              caption = paste("Top", n, "markers per cluster")) %>% show
    }
    markers
})
