#' Plot Known Markers
#'
#' @rdname plot_known_markers
#'
#' @param y [known_markers_detected()] [tibble] grouped by cluster.
#' @param markdown Print Markdown headers.
#'
#' @export
setMethod(
    "plot_known_markers",
    signature(x = "seurat", y = "grouped_df"),
    function(x, y, markdown = TRUE) {
        if (nrow(y) == 0L) {
            return(NULL)
        }
        cell_types <- y %>%
            pull("cell_type") %>%
            unique
        pblapply(seq_along(cell_types), function(a) {
            cell_type <- cell_types[[a]]
            if (isTRUE(markdown)) {
                writeLines(c(
                    "",
                    "",
                    paste("###", cell_type),
                    ""))
            }
            symbols <- y %>%
                filter(.data[["cell_type"]] == !!cell_type) %>%
                pull("symbol") %>%
                unique %>%
                sort
            if (!is.null(symbols)) {
                plot_clusters(x, symbols)
            } else {
                NULL
            }
        }) %>%
            invisible
    })
