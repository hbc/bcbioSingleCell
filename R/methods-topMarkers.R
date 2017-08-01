#' Top Markers
#'
#' @rdname topMarkers
#' @name topMarkers
#'
#' @param n Number of genes per cluster.
#' @param show Show [kable].
#'
#' @return [tibble].
#' @export
NULL



# Constructors ====
# Currently supports [Seurat::FindAllMarkers()]
.groupMarkers <- function(df) {
    if (!is.data.frame(df)) stop("data.frame required")
    # Rename `gene` to `symbol` if necessary
    if ("gene" %in% colnames(df)) {
        df <- rename(df, symbol = .data[["gene"]])
    }
    df %>%
        remove_rownames %>%
        as("tibble") %>%
        camel %>%
        tidy_select(c("cluster", "symbol"), everything()) %>%
        group_by(.data[["cluster"]])
}



.topMarkers <- function(object, n = 4L, show = FALSE) {
    markers <- .groupMarkers(object) %>%
        top_n(n = n, wt = .data[["avgDiff"]])
    if (isTRUE(show)) {
        kable(markers,
              caption = paste("Top", n, "markers per cluster")) %>% show
    }
    markers
}



# Methods ====
#' @rdname topMarkers
#' @export
setMethod("topMarkers", "data.frame", .topMarkers)
