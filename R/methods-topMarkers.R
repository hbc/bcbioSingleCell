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
.groupMarkers <- function(df) {
    df <- camel(df)

    # [Seurat::FindAllMarkers()] output
    if (identical(colnames(df),
        c("pVal", "avgDiff", "pct1", "pct2", "cluster", "gene"))) {
        # Rename the gene symbol and P value columns
        message("Fixing columns in Seurat marker data.frame")
        df <- df %>%
            rename(pvalue = .data[["pVal"]],
                   symbol = .data[["gene"]])
    }

    # Ensure that required columns are present
    requiredCols <- c("cluster", "ensgene", "symbol", "pvalue")
    if (!all(c("cluster", "ensgene", "symbol", "pvalue") %in% colnames(df))) {
        stop(paste("Marker data.frame must contain:", toString(requiredCols)))
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
