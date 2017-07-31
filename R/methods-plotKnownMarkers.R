#' Plot Known Markers
#'
#' @rdname plotKnownMarkers
#' @author Michael Steinbaugh
#'
#' @param y [knownMarkersDetected()] [tibble] grouped by cluster.
#' @param markdown Print Markdown headers.
#'
#' @return [writeLines()].
#' @export
setMethod(
    "plotKnownMarkers",
    signature(x = "seurat", y = "grouped_df"),
    function(x, y, markdown = TRUE) {
        if (nrow(y) == 0L) {
            return(NULL)
        }
        cellTypes <- y %>%
            pull("cellType") %>%
            unique
        pblapply(seq_along(cellTypes), function(a) {
            cellType <- cellTypes[[a]]
            if (isTRUE(markdown)) {
                writeLines(c(
                    "",
                    "",
                    paste("###", cellType),
                    ""))
            }
            symbols <- y %>%
                filter(.data[["cellType"]] == !!cellType) %>%
                pull("symbol") %>%
                unique %>%
                sort
            if (!is.null(symbols)) {
                plotClusters(x, symbols)
            } else {
                NULL
            }
        }) %>%
            invisible
    })
