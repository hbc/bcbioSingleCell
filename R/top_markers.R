#' Top marker genes
#'
#' @param markers Markers.
#' @param n Number of genes per cluster.
#'
#' @return [data.frame].
#' @export
top_markers <- function(markers, n = 4L) {
    markers %>%
        group_by(.data[["cluster"]]) %>%
        top_n(n = n, wt = .data[["avg_diff"]])
}
