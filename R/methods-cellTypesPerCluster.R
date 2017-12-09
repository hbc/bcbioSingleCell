#' Cell Types per Cluster
#'
#' @rdname cellTypesPerCluster
#' @name cellTypesPerCluster
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param min Minimum number of marker genes per cluster.
#' @param max Maximum number of marker genes per cluster.
#'
#' @return [tibble] grouped by cluster, containing the count (`n`) of
#'   significant known makers per cell type.
#'
#' @examples
#' \dontrun{
#' knownMarkers <- knownMarkersDetected(
#'     all = seuratAllMarkers,
#'     known = cellTypeMarkers)
#' cellTypesPerCluster(knownMarkers)
#' }
NULL



# Constructors =================================================================
#' @importFrom dplyr desc everything group_by select ungroup
#' @importFrom rlang !!! quos
.cellTypesPerCluster <- function(
    object,
    min = 2,
    max = Inf) {
    if (attr(object, "vars") != "cell") {
        stop("Markers tibble should be grouped by cell", call. = FALSE)
    }
    requiredCols <- c(
        "avgLogFC",  # Seurat v2.1
        "cell",      # bcbio
        "cluster",   # Seurat
        "ensgene",   # bcbio
        "padj",      # Seurat v2.1
        "symbol"     # bcbio
    )
    if (!all(requiredCols %in% colnames(object))) {
        stop(paste(
            "Required columns:", toString(sort(requiredCols))
        ), call. = FALSE)
    }
    groupCols <- syms(c("cluster", "cell"))
    tbl <- object %>%
        ungroup() %>%
        # Use only positive markers for this approach
        filter(.data[["avgLogFC"]] > 0) %>%
        select(!!!groupCols, everything()) %>%
        group_by(!!!groupCols) %>%
        arrange(.data[["padj"]], .by_group = TRUE) %>%
        summarize(
            n = n(),
            # Genes are arranged by P value
            symbol = toString(.data[["symbol"]]),
            ensgene = toString(.data[["ensgene"]])
        ) %>%
        group_by(!!sym("cluster")) %>%
        arrange(desc(.data[["n"]]), .by_group = TRUE)
    # Apply minimum and maximum gene cutoffs
    if (is.numeric(min) & min > 1) {
        tbl <- filter(tbl, .data[["n"]] >= min)
    }
    if (is.numeric(max) & max > 1) {
        tbl <- filter(tbl, .data[["n"]] <= max)

    }
    tbl
}



# Methods =====
#' @rdname cellTypesPerCluster
#' @export
setMethod(
    "cellTypesPerCluster",
    signature("grouped_df"),
    .cellTypesPerCluster)
