#' Top markers
#'
#' @rdname top_markers
#'
#' @param object Primary object.
#' @param n Number of genes per cluster.
#'
#' @return [tibble].
#' @export
setMethod("top_markers", "data.frame", function(object, n = 4L) {
    # Currently this supports [Seurat::FindAllMarkers()] return
    object %>%
        remove_rownames %>%
        as("tibble") %>%
        snake %>%
        tidy_select(c("cluster", "gene"), everything()) %>%
        group_by(.data[["cluster"]]) %>%
        top_n(n = n, wt = .data[["avg_diff"]])
})
